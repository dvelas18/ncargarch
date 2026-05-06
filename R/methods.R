#' @export
print.summary.ncar_fit <- function(x, ...) {
  cat(x$model, "\n", sep = "")
  cat("Distribution: ", x$distribution, "\n", sep = "")
  cat("\nStep 1: LS estimates\n")
  print(x$step1, row.names = FALSE)
  cat("\nStep 2: GARCH(1,1) QML estimates\n")
  print(x$step2, row.names = FALSE)
  cat("\nStep 3: WLS estimates\n")
  print(x$step3, row.names = FALSE)
  cat("\nGARCH convergence code: ", x$convergence, "\n", sep = "")
  invisible(x)
}

#' @export
print.summary.ncar_sim <- function(x, ...) {
  cat("Summary of ncar_sim object\n")
  cat("\nParameters\n")
  print(x$parameters)
  cat("\nDistribution\n")
  print(x$distribution)
  cat("\nMoments\n")
  print(x$moments, row.names = FALSE)
  invisible(x)
}
