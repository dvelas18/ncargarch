#' Plot estimated volatilities from an `ncar_fit` object
#'
#' @param object An `ncar_fit` object.
#' @param include_x Logical. If `TRUE`, include the AR(1) observed-process
#'   volatility measure when available.
#' @param ... Additional arguments passed to `plot()`.
#'
#' @return Invisibly returns `object`.
#' @export
ncar_plot_volatility <- function(object, include_x = TRUE, ...) {
  if (!inherits(object, "ncar_fit")) stop("object must be an ncar_fit object.", call. = FALSE)
  x_axis <- if (!is.null(object$dates_effective)) object$dates_effective else seq_along(object$sigma_epsilon)
  plot(x_axis, object$sigma_epsilon, type = "l", xlab = "Time", ylab = "Volatility", main = "Estimated volatility", ...)
  if (include_x && !all(is.na(object$sigma_x))) {
    lines(x_axis, object$sigma_x, lty = 2)
    legend("topright", legend = c("epsilon volatility", "X volatility"), lty = c(1, 2), bty = "n")
  }
  invisible(object)
}

#' Plot diagnostic ACFs from an `ncar_fit` object
#'
#' @param object An `ncar_fit` object.
#' @param lag.max Maximum lag passed to `stats::acf()`.
#'
#' @return Invisibly returns `object`.
#' @export
ncar_plot_diagnostics <- function(object, lag.max = 40L) {
  if (!inherits(object, "ncar_fit")) stop("object must be an ncar_fit object.", call. = FALSE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(2, 2))
  stats::acf(object$residuals, lag.max = lag.max, main = "ACF of residuals")
  stats::acf(object$residuals^2, lag.max = lag.max, main = "ACF of squared residuals")
  stats::acf(object$std_residuals, lag.max = lag.max, main = "ACF of standardized residuals")
  stats::acf(object$std_residuals^2, lag.max = lag.max, main = "ACF of squared standardized residuals")
  invisible(object)
}

#' @export
plot.ncar_fit <- function(x,
                          which = c("series", "volatility", "standardized", "diagnostics"),
                          ...) {
  which <- match.arg(which, several.ok = TRUE)
  if ("diagnostics" %in% which && length(which) == 1L) {
    return(ncar_plot_diagnostics(x))
  }

  panels <- setdiff(which, "diagnostics")
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(length(panels), 1L))
  x_axis <- if (!is.null(x$dates_effective)) x$dates_effective else seq_along(x$residuals)

  for (w in panels) {
    if (w == "series") {
      plot(x_axis, x$step1$y, type = "l", xlab = "Time", ylab = expression(X[t]), main = "Observed series", ...)
    }
    if (w == "volatility") {
      plot(x_axis, x$sigma_epsilon, type = "l", xlab = "Time", ylab = "Volatility", main = "Estimated innovation volatility", ...)
      if (!all(is.na(x$sigma_x))) lines(x_axis, x$sigma_x, lty = 2)
    }
    if (w == "standardized") {
      plot(x_axis, x$std_residuals, type = "l", xlab = "Time", ylab = expression(hat(eta)[t]), main = "Standardized residuals", ...)
      abline(h = 0, lty = 2)
    }
  }
  if ("diagnostics" %in% which) ncar_plot_diagnostics(x)
  invisible(x)
}
