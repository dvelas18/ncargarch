#' Simulate a noncausal AR(p) process with causal GARCH(1,1) volatility
#'
#' Simulates
#'
#' \deqn{X_t = \phi_1 X_{t+1} + \cdots + \phi_p X_{t+p} + \epsilon_t,}
#'
#' with
#'
#' \deqn{\epsilon_t = \sigma_t \eta_t, \quad
#'       \sigma_t^2 = \omega + \alpha \epsilon_{t-1}^2 + \beta \sigma_{t-1}^2.}
#'
#' The AR recursion is solved backward. A finite `future_padding` is used to
#' approximate the infinite noncausal representation.
#'
#' @param n Sample size to return.
#' @param phi Numeric vector of noncausal AR coefficients.
#' @param omega GARCH intercept.
#' @param alpha ARCH coefficient.
#' @param beta GARCH coefficient.
#' @param dist Innovation distribution: `"norm"`, `"std"`, or `"sstd"`.
#' @param df Degrees of freedom for `"std"` or `"sstd"`.
#' @param skew Skewness parameter for `"sstd"`.
#' @param burnin Number of initial observations discarded after the backward AR
#'   recursion.
#' @param future_padding Number of additional future observations used in the
#'   backward recursion. If `NULL`, it is set equal to `burnin`.
#' @param shock Optional shock. Use `FALSE`/`NULL` for no shock, or a list with
#'   fields `time`, `value`, and `type`. `type` can be `"eta"` or `"epsilon"`.
#'   The `time` is relative to the returned sample, so `time = 1` shocks the
#'   first returned observation.
#' @param seed Optional random seed.
#' @param return_burnin Logical. If `TRUE`, return the full simulated path in
#'   addition to the trimmed sample.
#'
#' @return An object of class `ncar_sim`.
#' @export
ncar_simulate <- function(n,
                          phi,
                          omega,
                          alpha,
                          beta,
                          dist = c("norm", "std", "sstd"),
                          df = 10,
                          skew = 1,
                          burnin = 500L,
                          future_padding = NULL,
                          shock = FALSE,
                          seed = NULL,
                          return_burnin = FALSE) {
  dist <- .ncar_match_dist(dist)
  phi <- as.numeric(phi)
  p <- length(phi)

  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 5) {
    stop("n must be a single integer greater than or equal to 5.", call. = FALSE)
  }
  n <- as.integer(n)
  burnin <- as.integer(burnin)
  if (burnin < 0L) stop("burnin must be non-negative.", call. = FALSE)
  if (is.null(future_padding)) future_padding <- burnin
  future_padding <- as.integer(future_padding)
  if (future_padding < p) future_padding <- p

  ncar_check_params(
    phi = phi,
    omega = omega,
    alpha = alpha,
    beta = beta,
    require_stationary_garch = TRUE,
    stop_on_error = TRUE
  )

  if (!is.null(seed)) set.seed(seed)

  total <- n + burnin + future_padding
  eta <- .ncar_draw_innovations(total + p, dist = dist, df = df, skew = skew)

  shock_info <- NULL
  shock_abs <- NA_integer_
  shock_type <- "none"
  shock_value <- NA_real_

  if (!(is.logical(shock) && length(shock) == 1L && !shock) && !is.null(shock)) {
    if (!is.list(shock)) stop("shock must be FALSE, NULL, or a list.", call. = FALSE)
    if (!all(c("time", "value") %in% names(shock))) {
      stop("shock must contain at least fields 'time' and 'value'.", call. = FALSE)
    }
    shock_type <- if ("type" %in% names(shock)) shock$type else "eta"
    shock_type <- match.arg(shock_type, choices = c("eta", "epsilon"))
    shock_time <- as.integer(shock$time)
    if (length(shock_time) != 1L || !is.finite(shock_time) || shock_time < 1L || shock_time > n) {
      stop("shock$time must be an integer between 1 and n, relative to the returned sample.", call. = FALSE)
    }
    shock_value <- as.numeric(shock$value)[1L]
    if (!is.finite(shock_value)) stop("shock$value must be finite.", call. = FALSE)
    shock_abs <- burnin + shock_time
    if (shock_type == "eta") eta[shock_abs] <- shock_value
    shock_info <- list(time = shock_time, absolute_time = shock_abs, value = shock_value, type = shock_type)
  }

  sigma2 <- numeric(total + p)
  eps <- numeric(total + p)
  sigma2[1L] <- omega / (1 - alpha - beta)
  eps[1L] <- sqrt(sigma2[1L]) * eta[1L]
  if (shock_type == "epsilon" && shock_abs == 1L) {
    eps[1L] <- shock_value
    eta[1L] <- eps[1L] / sqrt(sigma2[1L])
  }

  if (total + p >= 2L) {
    for (t in 2:(total + p)) {
      sigma2[t] <- omega + alpha * eps[t - 1L]^2 + beta * sigma2[t - 1L]
      if (shock_type == "epsilon" && t == shock_abs) {
        eps[t] <- shock_value
        eta[t] <- eps[t] / sqrt(sigma2[t])
      } else {
        eps[t] <- sqrt(sigma2[t]) * eta[t]
      }
    }
  }

  x <- numeric(total + p)
  x[(total + 1L):(total + p)] <- 0
  for (t in total:1L) {
    x[t] <- sum(phi * x[t + seq_len(p)]) + eps[t]
  }

  keep <- (burnin + 1L):(burnin + n)
  var_x <- if (p == 1L) .ncar_var_x_ar1(phi, omega, alpha, beta, sigma2[keep]) else rep(NA_real_, n)

  out <- list(
    x = x[keep],
    epsilon = eps[keep],
    eta = eta[keep],
    sigma_epsilon = sqrt(sigma2[keep]),
    var_epsilon = sigma2[keep],
    sigma_x = sqrt(var_x),
    var_x = var_x,
    params = list(phi = phi, omega = omega, alpha = alpha, beta = beta),
    distribution = list(dist = dist, df = df, skew = skew),
    shock = shock_info,
    settings = list(n = n, p = p, burnin = burnin, future_padding = future_padding, seed = seed),
    index = seq_len(n)
  )

  if (return_burnin) {
    out$full <- list(
      x = x[seq_len(total)],
      epsilon = eps[seq_len(total)],
      eta = eta[seq_len(total)],
      sigma_epsilon = sqrt(sigma2[seq_len(total)]),
      var_epsilon = sigma2[seq_len(total)]
    )
  }

  class(out) <- "ncar_sim"
  out
}

#' Simulate a noncausal AR(1)-GARCH(1,1) process
#'
#' Convenience wrapper around [ncar_simulate()] for the AR(1) case.
#'
#' @inheritParams ncar_simulate
#' @param phi Scalar noncausal AR coefficient.
#'
#' @return An object of class `ncar_sim`.
#' @export
ncar_simulate_ar1 <- function(n,
                              phi,
                              omega,
                              alpha,
                              beta,
                              dist = c("norm", "std", "sstd"),
                              df = 10,
                              skew = 1,
                              burnin = 500L,
                              future_padding = NULL,
                              shock = FALSE,
                              seed = NULL,
                              return_burnin = FALSE) {
  ncar_simulate(
    n = n,
    phi = phi,
    omega = omega,
    alpha = alpha,
    beta = beta,
    dist = dist,
    df = df,
    skew = skew,
    burnin = burnin,
    future_padding = future_padding,
    shock = shock,
    seed = seed,
    return_burnin = return_burnin
  )
}

#' Simulate a noncausal AR(p)-GARCH(1,1) process
#'
#' Convenience wrapper around [ncar_simulate()] for `length(phi) >= 1`.
#'
#' @inheritParams ncar_simulate
#'
#' @return An object of class `ncar_sim`.
#' @export
ncar_simulate_arp <- function(n,
                              phi,
                              omega,
                              alpha,
                              beta,
                              dist = c("norm", "std", "sstd"),
                              df = 10,
                              skew = 1,
                              burnin = 500L,
                              future_padding = NULL,
                              shock = FALSE,
                              seed = NULL,
                              return_burnin = FALSE) {
  ncar_simulate(
    n = n,
    phi = phi,
    omega = omega,
    alpha = alpha,
    beta = beta,
    dist = dist,
    df = df,
    skew = skew,
    burnin = burnin,
    future_padding = future_padding,
    shock = shock,
    seed = seed,
    return_burnin = return_burnin
  )
}

#' @export
print.ncar_sim <- function(x, ...) {
  cat("Noncausal AR(", x$settings$p, ") with causal GARCH(1,1) simulation\n", sep = "")
  cat("n = ", x$settings$n, ", dist = ", x$distribution$dist, "\n", sep = "")
  cat("phi = ", paste(format(x$params$phi, digits = 5), collapse = ", "), "\n", sep = "")
  cat("omega = ", format(x$params$omega, digits = 5),
      ", alpha = ", format(x$params$alpha, digits = 5),
      ", beta = ", format(x$params$beta, digits = 5), "\n", sep = "")
  if (!is.null(x$shock)) {
    cat("shock: type = ", x$shock$type, ", time = ", x$shock$time,
        ", value = ", x$shock$value, "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.ncar_sim <- function(object, ...) {
  out <- list(
    parameters = object$params,
    distribution = object$distribution,
    settings = object$settings,
    shock = object$shock,
    moments = data.frame(
      variable = c("x", "epsilon", "eta", "sigma_epsilon"),
      mean = c(mean(object$x), mean(object$epsilon), mean(object$eta), mean(object$sigma_epsilon)),
      sd = c(stats::sd(object$x), stats::sd(object$epsilon), stats::sd(object$eta), stats::sd(object$sigma_epsilon)),
      min = c(min(object$x), min(object$epsilon), min(object$eta), min(object$sigma_epsilon)),
      max = c(max(object$x), max(object$epsilon), max(object$eta), max(object$sigma_epsilon))
    )
  )
  class(out) <- "summary.ncar_sim"
  out
}

#' @export
plot.ncar_sim <- function(x,
                          which = c("x", "epsilon_volatility", "x_volatility"),
                          ...) {
  which <- match.arg(which, several.ok = TRUE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(length(which), 1L))

  for (w in which) {
    if (w == "x") {
      plot(x$index, x$x, type = "l", xlab = "Time", ylab = expression(X[t]), main = "Simulated noncausal process", ...)
      if (!is.null(x$shock)) abline(v = x$shock$time, lty = 2)
    }
    if (w == "epsilon_volatility") {
      plot(x$index, x$sigma_epsilon, type = "l", xlab = "Time", ylab = expression(sigma[epsilon,t]), main = "Innovation volatility", ...)
      if (!is.null(x$shock)) abline(v = x$shock$time, lty = 2)
    }
    if (w == "x_volatility") {
      if (all(is.na(x$sigma_x))) {
        plot.new(); title("X_t volatility is available only for AR(1) with alpha + beta < 1")
      } else {
        plot(x$index, x$sigma_x, type = "l", xlab = "Time", ylab = expression(sigma[X,t]), main = "Observed-process volatility under innovation filtration", ...)
        if (!is.null(x$shock)) abline(v = x$shock$time, lty = 2)
      }
    }
  }
  invisible(x)
}
