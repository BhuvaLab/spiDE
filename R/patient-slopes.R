# The slopes engine's per-patient niche slopes.
#
# From the shared, condition-free fit of index type k (patient intercepts,
# covariates gamma, pooled niche slopes beta; profile dispersion psi), every
# patient s gets a one-step negative binomial estimate of ITS OWN niche slopes,
# holding gamma and psi at the shared fit and profiling out its intercept:
#   delta_s = I_s^-1 u_s,   u_s = sum_{i in s} l~_i r_i,
#   I_s = sum_{i in s} w_i l~_i l~_i',   beta_s = beta + delta_s,
# with w = mu / (1 + psi mu), r = (y - mu) / (1 + psi mu) at the shared fit and
# l~ the niche log-densities centred within the patient with weights w. This is
# one Fisher-scoring step from the pooled slope toward the patient's own NB
# optimum: it uses each cell's count through the NB likelihood.
#
# Two variances per slope (research/simplify, 2026-09-28/29):
#   v_model = diag(I_s^-1), which treats a patient's cells as independent; and
#   v_tile, a within-patient SPATIAL sandwich I_s^-1 (sum_t U_t U_t') I_s^-1 *
#     T / (T - 1) over square tiles t of side `tile` in the patient's section(s),
#     U_t = sum_{i in t} l~_i r_i. Residual expression is spatially
#     autocorrelated, so for cell-rich patients v_model understates the slope's
#     variance. Patients with fewer than `min.tiles` tiles get v_model x the
#     gene's median tile/model ratio (at least 1).
# The condition test uses v_tile. The pooled test uses v_model x f_sn, a single
# per-(patient, niche) factor f_sn = max(1, median over genes of v_tile /
# v_model) (.patientFactor()): the per-gene tile sandwich is noisy on sparse
# genes and can sit near zero, which broke the pooled test's null.
#
# Genes are independent given the shared fit, so the loop is blocked and
# parallelised over genes; the per-patient factor is computed on the whole gene
# family afterwards, never within a block.

.patientSlopes <- function(fit, des, Yk, L, tiles, min.cells = 10L, min.tiles = 5L,
                           offset = NULL, BPPARAM = BiocParallel::SerialParam()) {
  genes <- rownames(Yk)
  G <- length(genes); S <- des$npat
  nt <- length(des$tested_niche)
  jn <- des$niche_cols                       # pooled slope columns in W, tested first
  rows_of <- split(seq_along(des$patient), factor(des$patient, levels = seq_len(S)))
  # each patient's absorbed block (its intercept, and under depth =
  # "spatial_spline" its library-size spline) over its own cells
  Z_of <- .patientBlocks(des, rows_of)
  one <- function(gi) {
    b <- v <- vm <- matrix(NA_real_, length(gi), S * nt)
    for (a in seq_along(gi)) {
      g <- gi[a]
      al <- fit$alpha[g, ]; psi <- fit$psi[g]
      if (!all(is.finite(al)) || !is.finite(psi)) next
      wr <- .workingWR(as.numeric(Yk[g, ]), des$W, al, psi,
                       offset = if (is.null(offset)) 0 else offset[g, ])
      w <- wr$w; r <- wr$r
      beta <- al[jn]
      bb <- vv <- vt <- matrix(NA_real_, S, nt)
      for (s in seq_len(S)) {
        i <- rows_of[[s]]
        if (length(i) < min.cells) next
        ws <- w[i]; Ls <- L[i, , drop = FALSE]
        Lt <- .partialBlock(Ls, ws, Z_of[[s]])
        if (is.null(Lt)) next
        Is <- crossprod(Lt * sqrt(ws))
        us <- colSums(Lt * r[i])
        # a non-finite weight or residual (an overflowing mean) would make
        # eigen() stop and take the whole gene block with it: this patient drops
        if (!all(is.finite(Is)) || !all(is.finite(us))) next
        e <- eigen(Is, symmetric = TRUE)
        # a patient whose niche columns are (near) collinear within its cells
        # gives no usable slope for this gene: it drops out, not the gene
        if (!all(is.finite(e$values)) || e$values[length(e$values)] <= max(e$values[1], 1e-12) * 1e-8) next
        Iinv <- e$vectors %*% (t(e$vectors) / e$values)
        d <- Iinv %*% us
        bb[s, ] <- (beta + d)[seq_len(nt)]
        vv[s, ] <- diag(Iinv)[seq_len(nt)]
        tt <- tiles[i]; nT <- length(unique(tt))
        if (nT >= min.tiles) {
          Ut <- rowsum(Lt * r[i], tt)
          Ut <- sweep(Ut, 2, colMeans(Ut))
          Vs <- Iinv %*% crossprod(Ut) %*% Iinv * nT / (nT - 1)
          vt[s, ] <- diag(Vs)[seq_len(nt)]
        }
      }
      rat <- apply(vt / vv, 2, stats::median, na.rm = TRUE)
      rat[!is.finite(rat)] <- 1
      miss <- is.na(vt) & is.finite(vv)
      vt[miss] <- (vv * matrix(pmax(rat, 1), S, nt, byrow = TRUE))[miss]
      b[a, ] <- as.numeric(bb); vm[a, ] <- as.numeric(vv); v[a, ] <- as.numeric(vt)
    }
    list(idx = gi, b = b, v = v, vm = vm)
  }
  res <- .bpGenes(G, one, BPPARAM)
  B <- V <- VM <- matrix(NA_real_, G, S * nt)
  for (x in res) {
    B[x$idx, ] <- x$b; V[x$idx, ] <- x$v; VM[x$idx, ] <- x$vm
  }
  dn <- list(genes, des$patients, des$tested_niche)
  list(beta = array(B, c(G, S, nt), dimnames = dn),
       v_tile = array(V, c(G, S, nt), dimnames = dn),
       v_model = array(VM, c(G, S, nt), dimnames = dn),
       pooled = matrix(fit$alpha[, jn[seq_len(nt)]], G, nt, dimnames = list(genes, des$tested_niche)),
       ncells = stats::setNames(lengths(rows_of), des$patients))
}

# The per-(patient, niche) variance factor of the pooled test, over the whole
# gene family: f_sn = max(1, median_g v_tile / v_model).
.patientFactor <- function(v_tile, v_model) {
  r <- v_tile / v_model
  r[!is.finite(r)] <- NA
  f <- pmax(apply(r, c(2, 3), stats::median, na.rm = TRUE), 1)
  f[!is.finite(f)] <- 1
  f
}
