# The sandwich engine's inference: a patient-clustered CR2 sandwich with the
# Bell-McCaffrey Satterthwaite degrees of freedom (Bell & McCaffrey 2002;
# Pustejovsky & Tipton 2018) on the fixed-effects NB fit of one index type,
# clusters = patients. Checked against clubSandwich::vcovCR(type = "CR2",
# inverse_var = TRUE) and coef_test(test = "Satterthwaite") on the working
# weighted lm of the same fit (tests/testthat/test-sandwich.R): SE within 0.3%,
# df to two decimals.
#
# In the working-weights space of the converged IRLS,
#   w_i = mu_i / (1 + psi mu_i),  r_i = (y_i - mu_i) / (1 + psi mu_i),
#   X~ = the non-absorbed columns centred within patient with weights w,
#   B = (X~' W X~)^-1,  G_s = sum_{i in s} w_i x~_i x~_i',  h_s = sum_{i in s} x~_i r_i.
# With G_s = U S U' (rank r_s), R_s = S^1/2 U' factors the patient's weighted
# block, the nonzero eigenvalues of the patient's hat block H_ss are those of
# R_s B R_s' = V L V', and the CR2-adjusted patient score is
#   u_s = R_s' V D V' S^-1/2 U' h_s,   D = diag((1 - L)^-1/2),
# so V_CR2 = B (sum_s u_s u_s') B -- never an n_s x n_s matrix. For coefficient
# j (c = e_j) the Bell-McCaffrey df is tr(Gm)^2 / ||Gm||_F^2 with
#   Gm = diag(||D V' R_s B c||^2) - Q' B Q,   Q[, s] = R_s' V D V' R_s B c.
# Genes are independent, so this is blocked and parallelised over genes.

.sandwichCR2 <- function(fit, des, Yk, offset = NULL, BPPARAM = BiocParallel::SerialParam()) {
  genes <- rownames(Yk)
  Xd <- des$W[, des$dense, drop = FALSE]
  jt <- match(des$tested, des$dense)
  pid <- des$patient; S <- des$npat
  rows_of <- split(seq_along(pid), factor(pid, levels = seq_len(S)))
  Z_of <- lapply(seq_len(S), function(s) des$W[rows_of[[s]], des$block_cols[[s]], drop = FALSE])
  intercept_only <- all(lengths(des$block_cols) == 1L)
  one <- function(gi) {
    out <- vector("list", length(gi))
    for (a in seq_along(gi)) {
      g <- gi[a]
      al <- fit$alpha[g, ]; psi <- fit$psi[g]
      if (!all(is.finite(al)) || !is.finite(psi)) next
      wr <- .workingWR(as.numeric(Yk[g, ]), des$W, al, psi,
                       offset = if (is.null(offset)) 0 else offset[g, ])
      w <- wr$w; r <- wr$r
      if (intercept_only) {
        sw <- rowsum(w, pid, reorder = TRUE)[, 1]
        xbar <- rowsum(Xd * w, pid, reorder = TRUE) / sw
        Xt <- Xd - xbar[pid, , drop = FALSE]
      } else {
        # residualise on each patient's absorbed block (Frisch-Waugh-Lovell);
        # the blocks nest within patients, so CR2 on the residualised design
        # equals CR2 on the full one (checked against clubSandwich in
        # tests/testthat/test-depth.R)
        Xt <- Xd
        for (s in seq_len(S)) {
          i <- rows_of[[s]]
          if (!length(i)) next
          xs <- .partialBlock(Xd[i, , drop = FALSE], w[i], Z_of[[s]])
          if (is.null(xs)) { Xt <- NULL; break }
          Xt[i, ] <- xs
        }
        if (is.null(Xt)) next
      }
      if (!all(is.finite(w)) || !all(is.finite(r))) next   # an overflowing mean: the gene drops out
      Info <- crossprod(Xt * sqrt(w))
      B <- tryCatch(solve(Info), error = function(e) NULL)
      if (!is.null(B) && !all(is.finite(B))) B <- NULL
      if (is.null(B)) next           # singular information for this gene: it drops out
      d <- ncol(Xd)
      Us <- matrix(0, d, S)
      fac <- vector("list", S)
      for (s in seq_len(S)) {
        i <- rows_of[[s]]
        if (!length(i)) next
        Xs <- Xt[i, , drop = FALSE]
        Gs <- crossprod(Xs * sqrt(w[i]))
        hs <- colSums(Xs * r[i])
        e <- eigen(Gs, symmetric = TRUE)
        keep <- e$values > max(e$values) * 1e-10
        if (!any(keep)) next
        Uu <- e$vectors[, keep, drop = FALSE]; sv <- e$values[keep]
        Rs <- sqrt(sv) * t(Uu)
        Ms <- Rs %*% B %*% t(Rs)
        ev <- eigen((Ms + t(Ms)) / 2, symmetric = TRUE)
        lam <- pmin(ev$values, 1 - 1e-8)
        Dg <- 1 / sqrt(1 - lam)
        Vv <- ev$vectors
        qe <- (t(Uu) %*% hs) / sqrt(sv)
        Us[, s] <- t(Rs) %*% (Vv %*% (Dg * (t(Vv) %*% qe)))
        fac[[s]] <- list(Rs = Rs, Vv = Vv, Dg = Dg)
      }
      Vcr2 <- B %*% tcrossprod(Us) %*% B
      df <- vapply(jt, function(j) {
        Bc <- B[, j]
        pn <- numeric(S); Q <- matrix(0, d, S)
        for (s in seq_len(S)) {
          f <- fac[[s]]
          if (is.null(f)) next
          z <- t(f$Vv) %*% (f$Rs %*% Bc)
          pn[s] <- sum((f$Dg * z)^2)
          Q[, s] <- t(f$Rs) %*% (f$Vv %*% (f$Dg * z))
        }
        Gm <- diag(pn, S) - t(Q) %*% B %*% Q
        sum(diag(Gm))^2 / sum(Gm^2)
      }, numeric(1))
      out[[a]] <- data.frame(gene = genes[g], niche = des$tested_niche,
                             estimate = unname(al[des$tested]),
                             se = sqrt(pmax(diag(Vcr2)[jt], 0)), df = df,
                             se_model = sqrt(diag(B)[jt]), psi = unname(psi),
                             stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
  }
  res <- .bpGenes(length(genes), one, BPPARAM)
  d <- do.call(rbind, res)
  if (is.null(d)) {
    d <- data.frame(gene = character(), niche = character(), estimate = numeric(), se = numeric(),
                    df = numeric(), se_model = numeric(), psi = numeric())
  }
  d$t <- d$estimate / d$se
  d$p <- 2 * stats::pt(-abs(d$t), d$df)
  d
}
