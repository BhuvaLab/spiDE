# Stage-2 machinery: patient pooling of per-core slopes, the DL tau2, the
# weighted moderated limma contrast, and the dropout diagnostics.

#' Precision-weighted pooling of per-core slopes within patients
#'
#' Cores share their patient's slope, so fixed-effect (1/v) pooling is
#' correct here (the between-patient component enters later, in tau2).
#' @noRd
.poolPatientSlopes <- function(B, V, sample2patient) {
  pats <- unique(unname(sample2patient[dimnames(B)[[3]]]))
  db <- dimnames(B)
  Bp <- Vp <- array(NA_real_, c(dim(B)[1:2], length(pats)),
                    dimnames = list(db[[1]], db[[2]], pats))
  for (p in pats) {
    ss <- dimnames(B)[[3]][sample2patient[dimnames(B)[[3]]] == p]
    b <- B[, , ss, drop = FALSE]; v <- V[, , ss, drop = FALSE]
    wt <- 1 / v
    wt[is.na(b) | !is.finite(wt)] <- NA
    swt <- apply(wt, 1:2, sum, na.rm = TRUE)
    num <- apply(b * wt, 1:2, sum, na.rm = TRUE)
    est <- num / swt
    est[swt == 0] <- NA_real_
    Bp[, , p] <- est
    vp <- 1 / swt
    vp[swt == 0] <- NA_real_
    Vp[, , p] <- vp
  }
  list(beta = Bp, var = Vp)
}
