#' Check parameters for a noncausal AR(p)-GARCH(1,1) model
#'
#' Checks the noncausal AR polynomial and the GARCH(1,1) restrictions.
#'
#' The noncausal AR(p) model is written as
#' 
#' \deqn{X_t = \phi_1 X_{t+1} + \cdots + \phi_p X_{t+p} + \epsilon_t.}
#'
#' The AR polynomial is \eqn{1 - \phi_1 z - \cdots - \phi_p z^p}. The
#' stationary noncausal solution requires roots outside the unit disk.
#'
#' @param phi Numeric scalar or vector of noncausal AR coefficients.
#' @param omega GARCH intercept. Must be strictly positive.
#' @param alpha ARCH coefficient. Must be non-negative.
#' @param beta GARCH coefficient. Must be non-negative.
#' @param require_stationary_garch Logical. If `TRUE`, require
#'   `alpha + beta < 1`. If `FALSE`, only warn when the condition fails.
#' @param tol Numerical tolerance.
#' @param stop_on_error Logical. If `TRUE`, stop when restrictions fail.
#'
#' @return A list with fields `ok`, `messages`, `ar_roots`, and `root_modulus`.
#' @export
ncar_check_params <- function(phi,
                              omega,
                              alpha,
                              beta,
                              require_stationary_garch = FALSE,
                              tol = 1e-8,
                              stop_on_error = TRUE) {
  messages <- character(0)

  phi <- as.numeric(phi)
  if (!length(phi) || any(!is.finite(phi))) {
    messages <- c(messages, "phi must be a finite numeric scalar or vector.")
  }

  if (!is.numeric(omega) || length(omega) != 1L || !is.finite(omega) || omega <= 0) {
    messages <- c(messages, "omega must be a single strictly positive finite number.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0) {
    messages <- c(messages, "alpha must be a single non-negative finite number.")
  }
  if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta) || beta < 0) {
    messages <- c(messages, "beta must be a single non-negative finite number.")
  }

  ar_roots <- numeric(0)
  root_modulus <- numeric(0)
  if (length(messages) == 0L) {
    ar_roots <- polyroot(c(1, -phi))
    root_modulus <- Mod(ar_roots)
    if (any(root_modulus <= 1 + tol)) {
      messages <- c(
        messages,
        "The AR polynomial 1 - phi_1 z - ... - phi_p z^p has at least one root inside or on the unit disk."
      )
    }

    if (alpha + beta >= 1) {
      msg <- "alpha + beta >= 1. The GARCH unconditional variance omega/(1-alpha-beta) is not finite."
      if (require_stationary_garch) {
        messages <- c(messages, msg)
      } else {
        warning(msg, call. = FALSE)
      }
    }
  }

  ok <- length(messages) == 0L
  if (!ok && stop_on_error) {
    stop(paste(messages, collapse = "\n"), call. = FALSE)
  }

  list(
    ok = ok,
    messages = messages,
    ar_roots = ar_roots,
    root_modulus = root_modulus,
    garch_persistence = alpha + beta
  )
}

.ncar_require_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this function. Install it with install.packages('", pkg, "').", call. = FALSE)
  }
}

.ncar_match_dist <- function(dist) {
  match.arg(dist, choices = c("norm", "std", "sstd"))
}

.ncar_draw_innovations <- function(n, dist = "norm", df = 10, skew = 1) {
  dist <- .ncar_match_dist(dist)
  if (dist == "norm") {
    return(stats::rnorm(n))
  }
  if (dist == "std") {
    if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 2) {
      stop("df must be larger than 2 for standardized Student-t innovations.", call. = FALSE)
    }
    if (requireNamespace("rugarch", quietly = TRUE)) {
      return(rugarch::rdist("std", n = n, mu = 0, sigma = 1, shape = df))
    }
    return(sqrt((df - 2) / df) * stats::rt(n, df = df))
  }
  if (dist == "sstd") {
    .ncar_require_namespace("rugarch")
    if (!is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 2) {
      stop("df must be larger than 2 for skew Student-t innovations.", call. = FALSE)
    }
    if (!is.numeric(skew) || length(skew) != 1L || !is.finite(skew) || skew <= 0) {
      stop("skew must be a strictly positive finite number for skew Student-t innovations.", call. = FALSE)
    }
    return(rugarch::rdist("sstd", n = n, mu = 0, sigma = 1, skew = skew, shape = df))
  }
}

.ncar_garch11_filter <- function(eps,
                                 omega,
                                 alpha,
                                 beta,
                                 initial = c("unconditional", "sample"),
                                 sigma2_0 = NULL,
                                 floor = 1e-12) {
  initial <- match.arg(initial)
  eps <- as.numeric(eps)
  n <- length(eps)
  if (n < 1L) stop("eps must contain at least one observation.", call. = FALSE)

  sigma2 <- numeric(n)
  if (!is.null(sigma2_0)) {
    sigma2[1L] <- as.numeric(sigma2_0)[1L]
  } else if (initial == "unconditional" && alpha + beta < 1) {
    sigma2[1L] <- omega / (1 - alpha - beta)
  } else {
    sigma2[1L] <- stats::var(eps)
    if (!is.finite(sigma2[1L]) || sigma2[1L] <= 0) sigma2[1L] <- mean(eps^2)
  }
  sigma2[1L] <- max(sigma2[1L], floor)

  if (n >= 2L) {
    for (t in 2:n) {
      sigma2[t] <- omega + alpha * eps[t - 1L]^2 + beta * sigma2[t - 1L]
      sigma2[t] <- max(sigma2[t], floor)
    }
  }
  sigma2
}

.ncar_var_x_ar1 <- function(phi, omega, alpha, beta, sigma2) {
  if (length(phi) != 1L) return(rep(NA_real_, length(sigma2)))
  if (!is.finite(alpha + beta) || alpha + beta >= 1) return(rep(NA_real_, length(sigma2)))
  den1 <- 1 - phi^2
  den2 <- 1 - phi^2 * (alpha + beta)
  if (den1 <= 0 || den2 <= 0) return(rep(NA_real_, length(sigma2)))
  sigma_bar <- omega / (1 - alpha - beta)
  out <- sigma_bar / den1 + (sigma2 - sigma_bar) / den2
  out[out < 0] <- NA_real_
  out
}

.ncar_safe_solve <- function(a, b = NULL, ridge = 1e-10) {
  out <- try(if (is.null(b)) solve(a) else solve(a, b), silent = TRUE)
  if (!inherits(out, "try-error") && all(is.finite(out))) return(out)
  a2 <- a + diag(ridge, nrow(a))
  if (is.null(b)) solve(a2) else solve(a2, b)
}

.ncar_make_lag_matrix <- function(x, p) {
  x <- as.numeric(x)
  n_total <- length(x)
  if (p < 1L || p >= n_total) stop("p must satisfy 1 <= p < length(x).", call. = FALSE)
  n_eff <- n_total - p
  y <- x[seq_len(n_eff)]
  z <- matrix(NA_real_, nrow = n_eff, ncol = p)
  for (j in seq_len(p)) {
    z[, j] <- x[(1L + j):(n_eff + j)]
  }
  colnames(z) <- paste0("lead", seq_len(p))
  list(y = y, z = z, n_eff = n_eff)
}
