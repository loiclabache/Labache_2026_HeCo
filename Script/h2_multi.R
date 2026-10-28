################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 13, 2025                                                               #
################################################################################


# Multivariate heritability estimator...........................................
################################################################################
h2_multi_R <- function(P, K, X, n_perm = 0L, seed = 1L) {
  P <- as.matrix(P)
  K <- as.matrix(K)
  X <- as.matrix(X)
  
  n <- nrow(X)
  if (!all(dim(P) == c(n, n)) || !all(dim(K) == c(n, n))) {
    stop("P, K, and X have incompatible dimensions.", call. = FALSE)
  }
  if (anyNA(P) || anyNA(K) || anyNA(X)) stop("P, K, and X cannot contain NA values.")
  if (!isTRUE(all.equal(P, t(P), tolerance = 1e-8))) stop("P must be symmetric.")
  if (!isTRUE(all.equal(K, t(K), tolerance = 1e-8))) stop("K must be symmetric.")
  if (qr(X)$rank < ncol(X)) stop("The covariate matrix X is rank deficient.")
  
  set.seed(seed)
  
  XtX_inv <- solve(crossprod(X))
  P0 <- diag(n) - X %*% XtX_inv %*% t(X)
  U <- svd(P0)$u[, seq_len(n - ncol(X)), drop = FALSE]
  
  Pp <- crossprod(U, P %*% U)
  Kp <- crossprod(U, K %*% U)
  n_resid <- nrow(Pp)
  
  kappa <- sum(diag(Kp %*% Kp)) / n_resid
  tau <- sum(diag(Kp)) / n_resid
  vK <- n_resid * (kappa - tau^2)
  if (!is.finite(vK) || abs(vK) < .Machine$double.eps) {
    stop("Degenerate kinship matrix after covariate projection.", call. = FALSE)
  }
  
  estimate_h2 <- function(K_current) {
    Qg <- K_current - tau * diag(n_resid)
    Qe <- kappa * diag(n_resid) - tau * K_current
    tg <- sum(diag(Qg %*% Pp)) / vK
    te <- sum(diag(Qe %*% Pp)) / vK
    estimate <- tg / (tg + te)
    pmin(1, pmax(0, estimate))
  }
  
  h2 <- estimate_h2(Kp)
  if (n_perm == 0L) return(list(h2 = h2, p_perm = NA_real_, h2_perm = NULL))
  
  h2_perm <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    K_perm <- if (i == 1L) Kp else {
      perm <- sample.int(n_resid)
      Kp[perm, perm, drop = FALSE]
    }
    h2_perm[i] <- estimate_h2(K_perm)
  }
  
  list(h2 = h2, p_perm = mean(h2_perm >= h2), h2_perm = h2_perm)
}

jackknife_family_se <- function(P, K, X, famid) {
  famid <- as.character(famid)
  if (length(famid) != nrow(X)) stop("famid must have one value per participant.")
  
  families <- unique(famid)
  h2_leave_one_family_out <- numeric(length(families))
  for (i in seq_along(families)) {
    keep <- famid != families[i]
    h2_leave_one_family_out[i] <- h2_multi_R(P[keep, keep, drop = FALSE],
                                             K[keep, keep, drop = FALSE],
                                             X[keep, , drop = FALSE],
                                             n_perm = 0L)$h2
  }
  
  m <- length(families)
  se <- sqrt((m - 1) / m * sum((h2_leave_one_family_out - mean(h2_leave_one_family_out))^2))
  list(se = se, estimates = h2_leave_one_family_out)
}
