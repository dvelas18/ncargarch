#' Extract model components from an `ncar_fit` object
#'
#' Returns a data frame containing the observed series aligned with the
#' noncausal regression, residuals, standardized residuals, estimated innovation
#' variance, estimated innovation volatility, and the AR(1) observed-process
#' volatility measure when available.
#'
#' @param object An `ncar_fit` object.
#' @param include_leads Logical. If `TRUE`, include the lead regressors
#'   `X_{t+1}, ..., X_{t+p}`.
#'
#' @return A data frame.
#' @export
ncar_components <- function(object, include_leads = TRUE) {
  if (!inherits(object, "ncar_fit")) stop("object must be an ncar_fit object.", call. = FALSE)
  n <- length(object$residuals)
  out <- data.frame(
    t = seq_len(n),
    x = object$step1$y,
    fitted_step1 = object$step1$fitted,
    fitted_step3 = object$step3$fitted,
    epsilon_hat = object$residuals,
    eta_hat = object$std_residuals,
    sigma_epsilon_hat = object$sigma_epsilon,
    var_epsilon_hat = object$var_epsilon,
    sigma_x_hat = object$sigma_x,
    var_x_hat = object$var_x,
    weights = object$step3$weights,
    row.names = NULL
  )

  if (!is.null(object$dates_effective)) {
    out <- cbind(date = object$dates_effective, out)
  }

  if (include_leads) {
    leads <- as.data.frame(object$step1$z)
    names(leads) <- paste0("x_lead", seq_len(object$p))
    insert_at <- if (!is.null(object$dates_effective)) 4L else 3L
    out <- cbind(out[, seq_len(insert_at - 1L), drop = FALSE], leads, out[, insert_at:ncol(out), drop = FALSE])
  }

  out
}

#' Export components from an `ncar_fit` object
#'
#' @param object An `ncar_fit` object.
#' @param file Output file path.
#' @param format Output format: `"csv"` or `"rds"`.
#' @param include_leads Logical passed to [ncar_components()].
#' @param ... Additional arguments passed to `utils::write.csv()` when
#'   `format = "csv"`.
#'
#' @return Invisibly returns the output file path.
#' @export
ncar_export <- function(object,
                        file,
                        format = c("csv", "rds"),
                        include_leads = TRUE,
                        ...) {
  format <- match.arg(format)
  comp <- ncar_components(object, include_leads = include_leads)
  if (format == "csv") {
    utils::write.csv(comp, file = file, row.names = FALSE, ...)
  } else {
    saveRDS(comp, file = file)
  }
  invisible(file)
}
