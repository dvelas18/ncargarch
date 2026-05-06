.ncar_K_inf <- function(r) {
  Pbar <- (-r + sqrt(r * r + 4 * r)) / 2
  (Pbar + r) / (Pbar + r + 1)
}

.ncar_filter_level <- function(x, K) {
  n <- length(x)
  mu <- numeric(n)
  mu[1L] <- x[1L]
  if (n >= 2L) {
    for (t in 2:n) {
      mu[t] <- mu[t - 1L] + K * (x[t] - mu[t - 1L])
    }
  }
  mu
}

.ncar_bic_from_summary <- function(sm) {
  p <- sm$df[1L]
  dfres <- sm$df[2L]
  n_eff <- p + dfres
  rss <- (sm$sigma^2) * dfres
  n_eff * log(rss / n_eff) + p * log(n_eff)
}

.ncar_as_summary_lm <- function(obj) {
  if (inherits(obj, "summary.lm")) return(obj)
  if (inherits(obj, "lm")) return(summary(obj))
  stop("ur.df@testreg is neither 'lm' nor 'summary.lm'.", call. = FALSE)
}

.ncar_make_calpha <- function(fit, tau_row, alpha) {
  cvals <- fit@cval[tau_row, c("1pct", "5pct", "10pct")]
  avals <- c(0.01, 0.05, 0.10)
  stats::approx(x = avals, y = as.numeric(cvals), xout = alpha, rule = 2)$y
}

.ncar_adf_bic <- function(z, deterministic, kmax, alpha) {
  best_bic <- Inf
  best_fit <- NULL
  best_k <- NA_integer_
  tau_row <- switch(deterministic,
    none = "tau1",
    drift = "tau2",
    trend = "tau3"
  )

  for (k in 0:kmax) {
    fit_k <- try(
      urca::ur.df(y = z, type = deterministic, lags = k, selectlags = "Fixed"),
      silent = TRUE
    )
    if (inherits(fit_k, "try-error")) next
    sm <- .ncar_as_summary_lm(fit_k@testreg)
    bval <- .ncar_bic_from_summary(sm)
    if (is.finite(bval) && bval < best_bic) {
      best_bic <- bval
      best_fit <- fit_k
      best_k <- k
    }
  }

  if (is.null(best_fit)) stop("All ADF regressions failed.", call. = FALSE)
  sm <- .ncar_as_summary_lm(best_fit@testreg)
  coefs <- sm$coefficients
  rn <- rownames(coefs)
  idx <- which(rn == "z.lag.1")
  if (length(idx) != 1L) {
    idx <- grep("\\bz(\\.|)lag(\\.|)1\\b|^z\\.lag\\.1$", rn, perl = TRUE)
    if (length(idx) == 0L) stop("Cannot locate z_{t-1} in the ADF regression.", call. = FALSE)
    idx <- idx[1L]
  }

  tau <- unname(coefs[idx, "t value"])
  c_alpha <- .ncar_make_calpha(best_fit, tau_row = tau_row, alpha = alpha)
  reject <- is.finite(tau) && is.finite(c_alpha) && tau <= c_alpha

  list(k = best_k, tau = tau, c_alpha = c_alpha, reject = reject, bic = best_bic)
}

#' One-sided local-level detrending with ADF-based smoothing selection
#'
#' Applies the one-sided local-level smoother
#'
#' \deqn{\mu_t = \mu_{t-1} + K_\infty(r)(p_t - \mu_{t-1}),}
#'
#' and selects the signal-to-noise ratio `r` using an ADF stationarity check on
#' the detrended series `z_t = p_t - mu_t`.
#'
#' @param x Numeric price or level series.
#' @param r_grid Positive grid of signal-to-noise ratios.
#' @param alpha ADF significance level. The `urca` critical values are available
#'   at 1%, 5%, and 10%; other values are linearly interpolated or extrapolated.
#' @param deterministic Deterministic component in the ADF regression:
#'   `"none"`, `"drift"`, or `"trend"`.
#' @param kmax Maximum ADF lag. If `NULL`, uses the Schwert rule.
#' @param selection Selection rule. `"min_reject"` chooses the smallest `r` that
#'   rejects the unit-root null. `"closest_left"` reproduces the earlier script
#'   rule: among rejecting values, choose the statistic closest to the left
#'   critical value.
#' @param return_scan_table Logical. If `TRUE`, include the full grid scan.
#'
#' @return An object of class `ncar_detrend`.
#' @export
ncar_detrend <- function(x,
                         r_grid = exp(seq(log(1e-6), log(10), length.out = 100L)),
                         alpha = 0.01,
                         deterministic = c("none", "drift", "trend"),
                         kmax = NULL,
                         selection = c("min_reject", "closest_left"),
                         return_scan_table = TRUE) {
  .ncar_require_namespace("urca")
  deterministic <- match.arg(deterministic)
  selection <- match.arg(selection)

  x <- as.numeric(x)
  n <- length(x)
  if (n < 25L) stop("x is too short; at least 25 observations are recommended.", call. = FALSE)
  if (any(!is.finite(x))) stop("x contains non-finite values.", call. = FALSE)
  if (any(r_grid <= 0) || any(!is.finite(r_grid))) stop("r_grid must contain strictly positive finite values.", call. = FALSE)
  r_grid <- sort(unique(as.numeric(r_grid)))

  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha)) stop("alpha must be a single finite number.", call. = FALSE)
  if (alpha < 0.01 || alpha > 0.10) {
    warning("alpha is outside [0.01, 0.10]; using linear extrapolation from urca critical values.", call. = FALSE)
  }

  if (is.null(kmax)) {
    kmax <- floor(12 * (n / 100)^(1 / 4))
    kmax <- max(0L, kmax)
  }
  kmax <- as.integer(kmax)

  scan <- lapply(r_grid, function(r) {
    K <- .ncar_K_inf(r)
    mu <- .ncar_filter_level(x, K)
    z <- x - mu
    adf <- .ncar_adf_bic(z, deterministic = deterministic, kmax = kmax, alpha = alpha)
    distance_left <- if (adf$reject) adf$c_alpha - adf$tau else NA_real_
    data.frame(
      r = r,
      K = K,
      half_life = log(0.5) / log(1 - K),
      k = adf$k,
      tau = adf$tau,
      c_alpha = adf$c_alpha,
      reject = adf$reject,
      distance_left = distance_left,
      bic = adf$bic
    )
  })
  scan_df <- do.call(rbind, scan)
  rownames(scan_df) <- NULL

  rejecting <- scan_df[scan_df$reject, , drop = FALSE]
  if (nrow(rejecting) == 0L) {
    choice_row <- scan_df[which.max(scan_df$r), , drop = FALSE]
    status <- "no_rejection_on_grid"
    note <- paste0("No r in the grid yields ADF rejection at ", 100 * alpha, "%. Returned max(r_grid).")
  } else if (selection == "min_reject") {
    choice_row <- rejecting[order(rejecting$r), , drop = FALSE][1L, , drop = FALSE]
    status <- "selected_minimum_rejecting_r"
    note <- paste0("Chose the smallest r whose ADF statistic rejects at ", 100 * alpha, "%.")
  } else {
    choice_row <- rejecting[order(rejecting$distance_left, rejecting$r), , drop = FALSE][1L, , drop = FALSE]
    status <- "selected_closest_left"
    note <- paste0("Chose the rejecting r whose ADF statistic is closest to the left critical value at ", 100 * alpha, "%.")
  }

  r_star <- choice_row$r
  K_star <- choice_row$K
  trend <- .ncar_filter_level(x, K_star)
  detrended <- x - trend
  adf_star <- .ncar_adf_bic(detrended, deterministic = deterministic, kmax = kmax, alpha = alpha)

  out <- list(
    r_star = r_star,
    K_star = K_star,
    half_life_star = log(0.5) / log(1 - K_star),
    trend = trend,
    detrended = detrended,
    mu = trend,
    z = detrended,
    adf = list(
      k_bic = adf_star$k,
      tau = adf_star$tau,
      c_alpha = adf_star$c_alpha,
      reject = adf_star$reject,
      alpha = alpha,
      deterministic = deterministic
    ),
    grid_scan = if (return_scan_table) scan_df else NULL,
    settings = list(alpha = alpha, deterministic = deterministic, kmax = kmax, n = n, selection = selection),
    status = status,
    note = note
  )
  class(out) <- "ncar_detrend"
  out
}

#' @export
print.ncar_detrend <- function(x, ...) {
  cat("One-sided local-level detrending with ADF selection\n")
  cat("Selected r*: ", format(x$r_star, digits = 6), "\n", sep = "")
  cat("K_infty(r*): ", format(x$K_star, digits = 6),
      "   Half-life: ", format(x$half_life_star, digits = 6), " periods\n", sep = "")
  cat("ADF (", x$adf$deterministic, ", alpha = ", x$adf$alpha, "): ",
      "k_BIC = ", x$adf$k_bic,
      ", tau = ", format(x$adf$tau, digits = 5),
      ", c_alpha = ", format(x$adf$c_alpha, digits = 5),
      ", reject = ", x$adf$reject, "\n", sep = "")
  cat("Status: ", x$status, "\n", sep = "")
  if (!is.null(x$note)) cat("Note: ", x$note, "\n", sep = "")
  invisible(x)
}

#' @export
plot.ncar_detrend <- function(x, ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(2, 1))
  plot(x$trend + x$detrended, type = "l", xlab = "Time", ylab = "Level", main = "Original series and one-sided trend", ...)
  lines(x$trend, lwd = 2)
  plot(x$detrended, type = "l", xlab = "Time", ylab = "Detrended", main = "Detrended series", ...)
  invisible(x)
}
