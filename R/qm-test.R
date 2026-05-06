.ncar_compute_rhat_vec <- function(eta_hat, m) {
  n <- length(eta_hat)
  if (m >= n) stop("m must be strictly smaller than the length of eta_hat.", call. = FALSE)
  if (m < 1L) stop("m must be at least 1.", call. = FALSE)
  r_hat <- numeric(m)
  for (h in seq_len(m)) {
    r_hat[h] <- sum(eta_hat[(h + 1L):n] * eta_hat[1L:(n - h)]) / n
  }
  names(r_hat) <- paste0("r_", seq_len(m))
  r_hat
}

.ncar_compute_dsigma2_dphi_ar1 <- function(X, eps, alpha, beta) {
  n <- length(eps)
  if (length(X) != n + 1L) stop("Length of X must be length(eps) + 1.", call. = FALSE)
  d <- numeric(n)
  d[1L] <- 0
  if (n >= 2L) {
    for (t in 2:n) {
      d[t] <- -2 * alpha * eps[t - 1L] * X[t] + beta * d[t - 1L]
    }
  }
  d
}

.ncar_compute_Phi <- function(X, eps, sigma, eta, alpha, beta, phi) {
  n_eps <- length(eps)
  if (length(sigma) != n_eps || length(eta) != n_eps) {
    stop("eps, sigma, and eta must all have length T-p.", call. = FALSE)
  }
  if (length(X) != n_eps + 1L) {
    stop("Length of X must be length(eps) + 1 for the AR(1) Qm statistic.", call. = FALSE)
  }
  d <- .ncar_compute_dsigma2_dphi_ar1(X, eps, alpha, beta)
  mean(-eta * X[2:(n_eps + 1L)] * d / (sigma^3))
}

.ncar_compute_iota_vec <- function(X, eps, sigma, eta, nu, phi, m) {
  n_eps <- length(eps)
  if (length(sigma) != n_eps || length(eta) != n_eps || length(nu) != n_eps) {
    stop("eps, sigma, eta, and nu must all have length T-1.", call. = FALSE)
  }
  if (length(X) != n_eps + 1L) {
    stop("Length of X must be length(eps) + 1.", call. = FALSE)
  }
  if (m < 1L || m >= n_eps) stop("m must satisfy 1 <= m < length(eps).", call. = FALSE)

  mu_sigma2 <- mean(X[2:(n_eps + 1L)]^2 / (sigma^2))
  iota <- numeric(m)

  for (i in seq_len(m)) {
    sum_h <- 0
    for (h in seq_len(n_eps - 1L)) {
      t_lower <- max(1L, 1L + i - h)
      t_upper <- n_eps - h
      if (t_upper < t_lower) next
      idx <- t_lower:t_upper
      term <- nu[idx] * eps[idx] * sigma[idx + h] * eta[idx + h - i]
      sum_h <- sum_h + (phi^(h - 1L)) * mean(term)
    }
    iota[i] <- sum_h / mu_sigma2
  }
  names(iota) <- paste0("iota_", seq_len(m))
  iota
}

.ncar_compute_z_vec <- function(X, eps, sigma, eta, alpha, beta, phi, m) {
  n <- length(eps)
  if (length(X) != n + 1L) stop("Length of X must be length(eps) + 1.", call. = FALSE)
  if (length(sigma) != n || length(eta) != n) stop("eps, sigma, and eta must all have length T-1.", call. = FALSE)
  if (m < 1L || m >= n) stop("m must satisfy 1 <= m < length(eps).", call. = FALSE)

  d <- .ncar_compute_dsigma2_dphi_ar1(X, eps, alpha, beta)
  sigma2 <- sigma^2
  z <- numeric(m)

  for (i in seq_len(m)) {
    t_idx <- (i + 1L):n
    e1 <- mean(sigma[t_idx] / sigma[t_idx - i])
    term2 <- eta[t_idx - i] * eta[t_idx] *
      (d[t_idx] / sigma2[t_idx] + d[t_idx - i] / sigma2[t_idx - i])
    e2 <- mean(term2)

    z[i] <- -(phi^(i - 1L)) * e1 - 0.5 * e2
  }
  names(z) <- paste0("z_", seq_len(m))
  z
}

.ncar_compute_sigma2_3wls <- function(X, eps, sigma2, Phi) {
  n <- length(eps)
  if (length(sigma2) != n) stop("sigma2 must have length T-1.", call. = FALSE)
  if (length(X) != n + 1L) stop("Length of X must be length(eps) + 1.", call. = FALSE)

  mu_sigma2 <- mean(X[2:(n + 1L)]^2 / sigma2)
  EX2 <- mean(X^2)
  nu <- 1 / sigma2 + Phi / EX2
  Y <- nu * eps * X[2:(n + 1L)]
  var_Y <- stats::var(Y)
  sigma2_3WLS <- var_Y / (mu_sigma2^2)

  list(
    sigma2_3WLS = sigma2_3WLS,
    mu_sigma2 = mu_sigma2,
    var_Y = var_Y,
    nu = nu
  )
}

.ncar_compute_Qm <- function(r_hat_m, z_m, iota_vec, sigma2_3WLS, n) {
  m <- length(r_hat_m)
  if (length(z_m) != m || length(iota_vec) != m) {
    stop("r_hat_m, z_m, and iota_vec must have the same length.", call. = FALSE)
  }

  Z <- cbind(z_m, iota_vec)
  zz <- sum(z_m^2)
  ii <- sum(iota_vec^2)
  zi <- sum(z_m * iota_vec)

  M <- matrix(c(
    zz, 1 + zi,
    1 + zi, ii - sigma2_3WLS
  ), nrow = 2L, byrow = TRUE)

  if (!is.finite(det(M)) || abs(det(M)) < 1e-8) {
    stop("Matrix M is numerically singular; cannot compute Qm.", call. = FALSE)
  }

  S <- diag(m) - Z %*% solve(M) %*% t(Z)
  Q_m <- as.numeric(n * t(r_hat_m) %*% S %*% r_hat_m)
  p_val <- stats::pchisq(Q_m, df = m, lower.tail = FALSE)

  list(
    statistic = Q_m,
    p_value = p_val,
    M = M,
    Z = Z,
    S = S,
    n = n,
    df = m
  )
}

.ncar_compute_qm_pieces <- function(X, eps, sigma2, eta, alpha, beta, phi, m) {
  n_eps <- length(eps)
  if (length(X) != n_eps + 1L) {
    stop("For the AR(1) Qm statistic, length(X) must be length(eps) + 1.", call. = FALSE)
  }
  if (length(sigma2) != n_eps || length(eta) != n_eps) {
    stop("sigma2 and eta must have the same length as eps.", call. = FALSE)
  }
  if (m < 1L || m >= n_eps) stop("m must satisfy 1 <= m < length(eps).", call. = FALSE)

  sigma <- sqrt(sigma2)
  Phi_hat <- .ncar_compute_Phi(X, eps, sigma, eta, alpha, beta, phi)
  res_sigma <- .ncar_compute_sigma2_3wls(X, eps, sigma2, Phi_hat)
  iota_vec <- .ncar_compute_iota_vec(X, eps, sigma, eta, res_sigma$nu, phi, m)
  z_vec <- .ncar_compute_z_vec(X, eps, sigma, eta, alpha, beta, phi, m)
  r_hat <- .ncar_compute_rhat_vec(eta, m)

  qres <- .ncar_compute_Qm(
    r_hat_m = r_hat,
    z_m = z_vec,
    iota_vec = iota_vec,
    sigma2_3WLS = res_sigma$sigma2_3WLS,
    n = n_eps
  )

  list(
    Phi = Phi_hat,
    sigma2_3WLS = res_sigma$sigma2_3WLS,
    mu_sigma2 = res_sigma$mu_sigma2,
    var_Y = res_sigma$var_Y,
    nu = res_sigma$nu,
    r_hat = r_hat,
    iota = iota_vec,
    z = z_vec,
    statistic = qres$statistic,
    p_value = qres$p_value,
    M = qres$M,
    Z = qres$Z,
    S = qres$S,
    n = qres$n,
    df = qres$df
  )
}

#' Portmanteau diagnostic test for standardized residuals
#'
#' Computes diagnostic portmanteau tests for standardized residuals. The default
#' `type = "qm"` implements the modified Qm statistic for the AR(1)-causal
#' GARCH(1,1) model. Benchmark Box-Pierce and Ljung-Box tests are also available.
#'
#' @param object An `ncar_fit` object.
#' @param m A scalar or vector of lag truncation values.
#' @param M If supplied and `m` is `NULL`, compute tests for all lags
#'   `1, ..., M`.
#' @param type Test type: `"qm"`, `"boxpierce"`, or `"ljungbox"`.
#'
#' @return An object of class `ncar_test` containing a table and, for `type="qm"`,
#'   the internal estimated pieces.
#' @export
ncar_test <- function(object,
                      m = NULL,
                      M = NULL,
                      type = c("qm", "boxpierce", "ljungbox")) {
  if (!inherits(object, "ncar_fit")) stop("object must be an ncar_fit object.", call. = FALSE)
  type <- match.arg(type)

  if (is.null(m)) {
    if (!is.null(M)) {
      m <- seq_len(as.integer(M))
    } else {
      m <- 10L
    }
  }
  m <- unique(as.integer(m))
  if (any(!is.finite(m)) || any(m < 1L)) stop("m must contain positive integers.", call. = FALSE)

  eta <- object$std_residuals
  n <- length(eta)
  if (any(m >= n)) stop("All m values must be smaller than the standardized residual length.", call. = FALSE)

  rows <- vector("list", length(m))
  pieces <- vector("list", length(m))

  for (i in seq_along(m)) {
    mi <- m[i]
    if (type == "qm") {
      if (object$p != 1L) {
        rows[[i]] <- data.frame(
          m = mi,
          statistic = NA_real_,
          df = mi,
          p_value = NA_real_,
          reject_5pct = NA,
          reject_10pct = NA,
          note = "The exact Qm implementation is currently restricted to AR(1)."
        )
        next
      }

      pars <- object$theta_used_for_standardization
      res <- try(
        .ncar_compute_qm_pieces(
          X = object$x,
          eps = object$residuals,
          sigma2 = object$var_epsilon,
          eta = object$std_residuals,
          alpha = pars$alpha,
          beta = pars$beta,
          phi = object$step3$phi[1L],
          m = mi
        ),
        silent = TRUE
      )

      if (inherits(res, "try-error")) {
        rows[[i]] <- data.frame(
          m = mi,
          statistic = NA_real_,
          df = mi,
          p_value = NA_real_,
          reject_5pct = NA,
          reject_10pct = NA,
          note = conditionMessage(attr(res, "condition"))
        )
      } else {
        rows[[i]] <- data.frame(
          m = mi,
          statistic = res$statistic,
          df = mi,
          p_value = res$p_value,
          reject_5pct = res$p_value < 0.05,
          reject_10pct = res$p_value < 0.10,
          note = ""
        )
        pieces[[i]] <- res
      }
    } else {
      bt_type <- if (type == "boxpierce") "Box-Pierce" else "Ljung-Box"
      bt <- stats::Box.test(eta, lag = mi, type = bt_type, fitdf = 0)
      rows[[i]] <- data.frame(
        m = mi,
        statistic = as.numeric(bt$statistic),
        df = mi,
        p_value = as.numeric(bt$p.value),
        reject_5pct = as.numeric(bt$p.value) < 0.05,
        reject_10pct = as.numeric(bt$p.value) < 0.10,
        note = "benchmark test"
      )
    }
  }

  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  names(pieces) <- paste0("m", m)

  out <- list(
    type = type,
    table = table,
    pieces = pieces,
    settings = list(m = m, n = n, p = object$p)
  )
  class(out) <- "ncar_test"
  out
}

#' @export
print.ncar_test <- function(x, ...) {
  cat("ncar portmanteau test: type = ", x$type, "\n", sep = "")
  print(x$table, row.names = FALSE)
  invisible(x)
}
