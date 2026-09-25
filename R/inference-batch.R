# Batched / accelerator-aware helpers for the inference stage: per-gene Gram
# matrices, their batched inverse sub-blocks, and the two independent memory
# budgets that bound them. Split out of inference.R, which keeps the
# statistical path (Wald covariance, within-gene combination, blocked driver).
#
# Everything here is shape-and-memory plumbing that must behave identically on
# a base R matrix and a torch tensor -- the two backends differ in ways that
# are easy to get wrong (see .subsetBatch() on advanced indexing), which is why
# they are collected in one place rather than inlined at their call sites.

#' Select rows of a matrix or a torch tensor
#'
#' Row-subsetting helper for the gene sub-batching in
#' \code{.waldCauchyBlock()}: base R \code{[} does not carry over to torch
#' tensors (which need \code{torch_index_select()}), and the working weights
#' may be either depending on the backend.
#'
#' @param x a matrix or torch tensor.
#' @param ii integer row positions to keep.
#' @return \code{x} restricted to rows \code{ii}, same type as the input.
#' @noRd
.rowsOf <- function(x, ii) {
  if (SpaNorm::is_torch_tensor(x)) {
    idx <- torch::torch_tensor(as.integer(ii), dtype = torch::torch_long(),
                               device = x$device)
    return(torch::torch_index_select(x, 1, idx))
  }
  x[ii, , drop = FALSE]
}

#' Subset a (block, p, p) array/tensor to (block, k, k) via a column index
#'
#' On a base R array, \code{varcov[, sel, sel]} already takes the
#' \code{sel x sel} cross-block (R's array indexing is orthogonal per
#' dimension). On a torch tensor this is NOT true -- tensor \code{[}-indexing
#' with two non-scalar index vectors pairs them elementwise (NumPy/PyTorch
#' "advanced indexing" semantics), so the two dimensions must be subset one
#' at a time via \code{torch_index_select()} instead.
#'
#' @param varcov a \code{(block, p, p)} array or torch tensor.
#' @param sel integer positions to keep along both of the last two dimensions.
#' @return a \code{(block, k, k)} array or torch tensor, \code{k = length(sel)}.
#' @noRd
.subsetBatch <- function(varcov, sel) {
  if (SpaNorm::is_torch_tensor(varcov)) {
    sel_t <- torch::torch_tensor(as.integer(sel), dtype = torch::torch_long(),
                                 device = varcov$device)
    varcov <- torch::torch_index_select(varcov, 2, sel_t)
    return(torch::torch_index_select(varcov, 3, sel_t))
  }
  varcov[, sel, sel, drop = FALSE]
}

#' Per-slice quadratic form w'Vw of a (block, p, p) array/tensor
#'
#' Needed for the patient-level contrast: the weighted-average response effect
#' across cell types is w'alpha, whose variance is w'Vw. V exists inside the
#' Wald block but only its diagonal is retained, so the quadratic form must be
#' formed here before V is discarded.
#'
#' @param varcov a \code{(block, p, p)} array or torch tensor.
#' @param w a length-\code{p} contrast vector.
#' @return a length-\code{block} numeric vector of quadratic forms.
#' @noRd
.batchQuad <- function(varcov, w) {
  if (SpaNorm::is_torch_tensor(varcov)) {
    wt <- torch::torch_tensor(as.numeric(w), dtype = varcov$dtype,
                              device = varcov$device)
    vw <- torch::torch_matmul(varcov, wt$unsqueeze(1)$unsqueeze(3))
    return(as.numeric(torch::torch_sum(vw$squeeze(3) * wt$unsqueeze(1), dim = 2)))
  }
  apply(varcov, 1, function(V) as.numeric(crossprod(w, V %*% w)))
}

#' Extract the per-slice diagonal of a (block, p, p) array/tensor
#'
#' @param varcov a \code{(block, p, p)} array or torch tensor.
#' @return a \code{block x p} matrix (or torch tensor) of diagonals.
#' @noRd
.batchDiag <- function(varcov) {
  if (SpaNorm::is_torch_tensor(varcov)) {
    return(torch::torch_diagonal(varcov, dim1 = 2, dim2 = 3))
  }
  b <- dim(varcov)[1]
  p <- dim(varcov)[2]
  idx <- cbind(rep(seq_len(b), p), rep(seq_len(p), each = b),
              rep(seq_len(p), each = b))
  matrix(varcov[idx], nrow = b, ncol = p)
}


#' Number of genes per inference block under a GPU memory budget
#'
#' Picks \code{block.size} for the GPU inference path so the per-block
#' genes x cells intermediates plus the batched covariance intermediates stay
#' within the accelerator's memory budget. Returns \code{NULL} (a single
#' block, the pre-existing behaviour) whenever the CPU path is in use -- so
#' the CPU path is never affected by GPU blocking, exactly mirroring
#' \code{SpaNorm::geneBlockCount()}'s own invariant.
#'
#' This sizes the \strong{NB math} only -- \code{Yb}, \code{mub}, \code{wtb},
#' the NB log-pmf internals and (mixed) the dispersion intermediates, all
#' \code{block x ncells} -- via \code{SPIDE_TENSOR_MULT_NB}. The covariance
#' side is \emph{not} sized here: \code{.waldCauchyBlock()} sub-batches the
#' \code{(batch, p, p)} Gram/inverse stack internally (see
#' \code{.covBatchSize()}), so it is bounded whatever block size is chosen,
#' and the two concerns stay independent. The multiplier is a rough,
#' deliberately generous estimate that needs empirical validation on real
#' hardware, exactly as noted for SpaNorm's own
#' \code{GPU_BLOCK_TENSOR_MULTIPLIER}.
#'
#' @param ng number of genes.
#' @param ncells number of cells.
#' @param p_eff the covariance dimension (\code{ncol(Wsub)} fixed effects,
#'   \code{ncol(W_full)} mixed); retained for signature stability and future
#'   use -- the covariance term is now bounded by \code{.covBatchSize()}.
#' @param backend the resolved backend.
#' @param gpu.mem.budget \code{NULL} (auto-detect via
#'   \code{SpaNorm::getGPUMemoryBudget()}) or a budget in bytes.
#' @return a single integer block size, or \code{NULL} for a single block.
#' @noRd
SPIDE_TENSOR_MULT_NB <- 12
SPIDE_TENSOR_MULT_COV <- 6
# fraction of the memory budget each of the two independently-bounded stages
# (NB math per gene block, covariance per sub-batch) may claim
SPIDE_BUDGET_FRACTION <- 0.5
# default cap on the covariance stack when there is no GPU budget to consult
# (the CPU path); override with options(spiDE.cov.mem.budget = <bytes>)
SPIDE_COV_MEM_BUDGET_CPU <- 2e9

.inferenceBlockSize <- function(ng, ncells, p_eff, backend,
                                gpu.mem.budget = NULL) {
  if (!(backend %in% c("gpu", "auto") && SpaNorm::checkGPU())) {
    return(NULL)
  }
  budget <- SpaNorm::getGPUMemoryBudget(gpu.mem.budget)
  if (!is.finite(budget)) {
    return(NULL)
  }
  bytes <- SpaNorm::gpuDtypeBytes()
  per_gene <- bytes * as.numeric(ncells) * SPIDE_TENSOR_MULT_NB
  b <- max(1L, as.integer(floor(budget * SPIDE_BUDGET_FRACTION / per_gene)))
  if (b >= ng) NULL else b
}

#' Number of genes per covariance sub-batch
#'
#' Bounds the \code{(batch, p, p)} Gram/inverse stack (and, on the GPU, the
#' \code{(batch, ncells, p)} weighted design feeding it) inside
#' \code{.waldCauchyBlock()}. This is the knob that makes wide mixed-effects
#' designs tractable: peak covariance memory is linear in \code{p} and scales
#' with this batch, so a design too wide to process all at once is handled by
#' shrinking the sub-batch rather than failing.
#'
#' Applies on \strong{both} backends. The CPU path needs it just as much as
#' the GPU path -- a single block of 13,348 genes at \code{p = 4906} would
#' otherwise try to allocate a \code{(13348, 4906, 4906)} array -- and unlike
#' \code{.inferenceBlockSize()} it therefore does not return NULL for CPU.
#'
#' @param ncells number of cells.
#' @param p the covariance dimension (\code{ncol} of the Gram design).
#' @param backend the resolved backend.
#' @param gpu.mem.budget \code{NULL} (auto-detect) or a budget in bytes; only
#'   consulted on the GPU path.
#' @param nworkers how many workers will evaluate this concurrently. The budget
#'   is a figure for the machine (or the device), but \code{.waldCauchyBlock()}
#'   runs inside \code{bplapply()} and each forked worker claims it
#'   independently: at 64 workers the CPU default of 2e9 is a 128 GB claim, in
#'   a stage that has already been OOM-killed once at 503 GB MaxRSS. Dividing
#'   here makes the documented budget the total it says it is. The GPU path
#'   takes the same division -- one device, several processes.
#' @return a single integer, genes per covariance sub-batch (at least 1).
#' @noRd
.covBatchSize <- function(ncells, p, backend, gpu.mem.budget = NULL,
                          nworkers = 1L) {
  gpu_active <- backend %in% c("gpu", "auto") && SpaNorm::checkGPU()
  budget <- if (gpu_active) {
    SpaNorm::getGPUMemoryBudget(gpu.mem.budget)
  } else {
    getOption("spiDE.cov.mem.budget", SPIDE_COV_MEM_BUDGET_CPU)
  }
  if (!is.finite(budget)) {
    budget <- SPIDE_COV_MEM_BUDGET_CPU
  }
  budget <- budget / max(1L, as.integer(nworkers))
  bytes <- if (gpu_active) SpaNorm::gpuDtypeBytes() else 8
  ncells <- as.numeric(ncells)
  p <- as.numeric(p)

  # GPU: the (batch, ncells, p) weighted design plus its transpose view, then
  # the (batch, p, p) Gram/Cholesky/inverse stack. CPU: the (batch, p, p)
  # stack only -- construction there is one p x p crossprod at a time.
  per_gene <- if (gpu_active) {
    bytes * (ncells * p * 2 + p^2 * SPIDE_TENSOR_MULT_COV)
  } else {
    bytes * p^2 * SPIDE_TENSOR_MULT_COV
  }
  max(1L, as.integer(floor(budget * SPIDE_BUDGET_FRACTION / per_gene)))
}
