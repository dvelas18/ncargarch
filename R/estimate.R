.ncar_lse <- function(x, p) {
  lmats <- .ncar_make_lag_matrix(x, p)
  y <- lmats$y
  z <- lmats$z
  fit <- stats::lm.fit(z, y)
  phi <- as.numeric(fit$coefficients)
  names(phi) <- paste0("phi", seq_len(p))
  resid <- as.numeric(y - z %*% phi)
  df_resid <- max(length(y) - p, 1L)
  s2 <- sum(resid^2) / df_resid
  xtx_inv <- .ncar_safe_solve(crossprod(z))
  vc <- s2 * xtx_inv
  se <- sqrt(pmax(diag(vc), 0))
  names(se) <- names(phi)

  list(
    phi = phi,
    residuals = resid,
    fitted = as.numeric(z %*% phi),
    se = se,
    vcov = vc,
    t = phi / se,
    y = y,
    z = z,
    n_eff = length(y),
    sigma2 = s2
  )
}

.ncar_wlse <- function(y, z, sigma2, p, se = c("basic", "none")) {
  se <- match.arg(se)
  w <- 1 / as.numeric(sigma2)
  if (length(w) != length(y)) stop("sigma2 length must match the effective sample length.", call. = FALSE)
  if (any(!is.finite(w)) || any(w <= 0)) stop("All WLS weights must be finite and positive.", call. = FALSE)

  zw <- z * sqrt(w)
  yw <- y * sqrt(w)
  a <- crossprod(zw)
  b <- crossprod(zw, yw)
  phi <- as.numeric(.ncar_safe_solve(a, b))
  names(phi) <- paste0("phi", seq_len(p))
  fitted <- as.numeric(z %*% phi)
  resid <- as.numeric(y - fitted)

  vc <- matrix(NA_real_, nrow = p, ncol = p)
  rownames(vc) <- colnames(vc) <- names(phi)
  se_out <- rep(NA_real_, p)
  names(se_out) <- names(phi)

  if (se == "basic") {
    n <- length(y)
    a_n <- a / n
    meat <- crossprod(z * (w * resid)) / n
    vc <- .ncar_safe_solve(a_n) %*% meat %*% .ncar_safe_solve(a_n) / n
    rownames(vc) <- colnames(vc) <- names(phi)
    se_out <- sqrt(pmax(diag(vc), 0))
  }

  list(
    phi = phi,
    residuals = resid,
    fitted = fitted,
    se = se_out,
    vcov = vc,
    t = phi / se_out,
    weights = w
  )
}

.ncar_garch_spec <- function(dist = "norm", fixed.pars = NULL) {
  .ncar_require_namespace("rugarch")
  dist <- .ncar_match_dist(dist)
  rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = dist,
    fixed.pars = fixed.pars
  )
}

.ncar_fit_garch11 <- function(resid,
                              dist = "norm",
                              solver = "hybrid",
                              solver.control = list(),
                              fit.control = list(),
                              ...) {
  .ncar_require_namespace("rugarch")
  spec <- .ncar_garch_spec(dist = dist)
  fit <- try(
    rugarch::ugarchfit(
      spec = spec,
      data = as.numeric(resid),
      solver = solver,
      solver.control = solver.control,
      fit.control = fit.control,
      ...
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    stop("rugarch::ugarchfit failed: ", conditionMessage(attr(fit, "condition")), call. = FALSE)
  }

  coef_all <- rugarch::coef(fit)
  sigma2 <- as.numeric(rugarch::sigma(fit)^2)

  matcoef <- NULL
  se <- rep(NA_real_, length(coef_all))
  names(se) <- names(coef_all)
  if (!is.null(fit@fit$matcoef)) {
    matcoef <- fit@fit$matcoef
    if ("Std. Error" %in% colnames(matcoef)) {
      idx <- intersect(names(coef_all), rownames(matcoef))
      se[idx] <- matcoef[idx, "Std. Error"]
    }
  }

  list(
    coef = coef_all,
    se = se,
    matcoef = matcoef,
    sigma2 = sigma2,
    sigma = sqrt(sigma2),
    fit = fit,
    convergence = fit@fit$convergence,
    message = fit@fit$message
  )
}

.ncar_extract_garch_pars <- function(garch_fit) {
  cf <- garch_fit$coef
  required <- c("omega", "alpha1", "beta1")
  if (!all(required %in% names(cf))) {
    stop("The fitted rugarch object does not contain omega, alpha1, and beta1.", call. = FALSE)
  }
  list(
    omega = as.numeric(cf["omega"]),
    alpha = as.numeric(cf["alpha1"]),
    beta = as.numeric(cf["beta1"]),
    coef = cf
  )
}

.ncar_step_table <- function(estimate, se, statistic_name = "t_stat") {
  out <- data.frame(
    parameter = names(estimate),
    estimate = as.numeric(estimate),
    std_error = as.numeric(se),
    statistic = as.numeric(estimate / se),
    row.names = NULL
  )
  names(out)[4L] <- statistic_name
  out
}

#' Three-step estimation of a noncausal AR(p)-causal GARCH(1,1) model
#'
#' Estimates the model
#'
#' \deqn{X_t = \phi_1 X_{t+1} + \cdots + \phi_p X_{t+p} + \epsilon_t,}
#'
#' \deqn{\epsilon_t = \sigma_t \eta_t, \quad
#'       \sigma_t^2 = \omega + \alpha \epsilon_{t-1}^2 + \beta \sigma_{t-1}^2.}
#'
#' The function implements the three-step procedure: first-step least squares
#' for the noncausal AR parameter, second-step GARCH(1,1) QML on first-step
#' residuals, and third-step weighted least squares using the estimated GARCH
#' variances as weights.
#'
#' @param x Numeric vector of observations.
#' @param p Noncausal AR order.
#' @param dist Conditional innovation distribution passed to `rugarch`:
#'   `"norm"`, `"std"`, or `"sstd"`.
#' @param demean Logical. If `TRUE`, subtract the sample mean before estimation.
#' @param dates Optional date/index vector with the same length as `x`.
#' @param se Standard-error option. `"basic"` uses practical LS/WLS standard
#'   errors. `"theory"` uses the AR(1) theorem-based third-step standard error
#'   where available and falls back to `"basic"` otherwise.
#' @param solver Solver passed to `rugarch::ugarchfit()`.
#' @param solver.control Solver-control list passed to `rugarch::ugarchfit()`.
#' @param fit.control Fit-control list passed to `rugarch::ugarchfit()`.
#' @param refit_garch_after_wls Logical. If `TRUE`, refit the GARCH model to the
#'   third-step residuals. The default `FALSE` keeps the strict three-step
#'   definition: the GARCH parameters are estimated in Step 2.
#' @param ... Additional arguments passed to `rugarch::ugarchfit()`.
#'
#' @return An object of class `ncar_fit`.
#' @export
ncar_fit <- function(x,
                     p = 1L,
                     dist = c("norm", "std", "sstd"),
                     demean = FALSE,
                     dates = NULL,
                     se = c("basic", "theory"),
                     solver = "hybrid",
                     solver.control = list(),
                     fit.control = list(),
                     refit_garch_after_wls = FALSE,
                     ...) {
  dist <- .ncar_match_dist(dist)
  se <- match.arg(se)
  x_original <- as.numeric(x)
  if (any(!is.finite(x_original))) stop("x contains non-finite values.", call. = FALSE)
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 1L) stop("p must be a positive integer.", call. = FALSE)
  p <- as.integer(p)
  if (length(x_original) <= p + 5L) stop("x is too short for the requested AR order.", call. = FALSE)
  if (!is.null(dates) && length(dates) != length(x_original)) stop("dates must have the same length as x.", call. = FALSE)

  x_mean <- if (demean) mean(x_original) else 0
  x_work <- x_original - x_mean

  step1 <- .ncar_lse(x_work, p)
  garch1 <- .ncar_fit_garch11(
    resid = step1$residuals,
    dist = dist,
    solver = solver,
    solver.control = solver.control,
    fit.control = fit.control,
    ...
  )
  pars <- .ncar_extract_garch_pars(garch1)

  step3 <- .ncar_wlse(
    y = step1$y,
    z = step1$z,
    sigma2 = garch1$sigma2,
    p = p,
    se = "basic"
  )

  garch_for_std <- garch1
  theta_for_std <- pars
  if (refit_garch_after_wls) {
    garch_refit <- .ncar_fit_garch11(
      resid = step3$residuals,
      dist = dist,
      solver = solver,
      solver.control = solver.control,
      fit.control = fit.control,
      ...
    )
    garch_for_std <- garch_refit
    theta_for_std <- .ncar_extract_garch_pars(garch_refit)
  }

  sigma2_step3 <- .ncar_garch11_filter(
    eps = step3$residuals,
    omega = theta_for_std$omega,
    alpha = theta_for_std$alpha,
    beta = theta_for_std$beta,
    initial = "unconditional"
  )
  eta_hat <- step3$residuals / sqrt(sigma2_step3)

  var_x <- if (p == 1L) {
    .ncar_var_x_ar1(step3$phi[1L], theta_for_std$omega, theta_for_std$alpha, theta_for_std$beta, sigma2_step3)
  } else {
    rep(NA_real_, length(step3$residuals))
  }

  theory <- NULL
  if (se == "theory" && p == 1L) {
    theory <- try(
      .ncar_compute_qm_pieces(
        X = x_work,
        eps = step3$residuals,
        sigma2 = sigma2_step3,
        eta = eta_hat,
        alpha = theta_for_std$alpha,
        beta = theta_for_std$beta,
        phi = step3$phi[1L],
        m = 1L
      ),
      silent = TRUE
    )
    if (!inherits(theory, "try-error") && is.finite(theory$sigma2_3WLS)) {
      step3$vcov[1L, 1L] <- theory$sigma2_3WLS / length(step3$residuals)
      step3$se[1L] <- sqrt(step3$vcov[1L, 1L])
      step3$t[1L] <- step3$phi[1L] / step3$se[1L]
    }
  }

  date_eff <- if (is.null(dates)) NULL else dates[seq_along(step3$residuals)]

  out <- list(
    x = x_work,
    x_original = x_original,
    x_mean = x_mean,
    p = p,
    dates = dates,
    dates_effective = date_eff,
    distribution = dist,
    step1 = step1,
    step2 = list(
      coef = garch1$coef,
      se = garch1$se,
      matcoef = garch1$matcoef,
      sigma2 = garch1$sigma2,
      sigma = garch1$sigma,
      fit = garch1$fit,
      convergence = garch1$convergence,
      message = garch1$message
    ),
    step3 = step3,
    residuals = step3$residuals,
    std_residuals = eta_hat,
    sigma_epsilon = sqrt(sigma2_step3),
    var_epsilon = sigma2_step3,
    sigma_x = sqrt(var_x),
    var_x = var_x,
    garch_refit_after_wls = if (refit_garch_after_wls) garch_for_std else NULL,
    theta_used_for_standardization = theta_for_std,
    settings = list(
      p = p,
      dist = dist,
      demean = demean,
      se = se,
      solver = solver,
      refit_garch_after_wls = refit_garch_after_wls
    ),
    theory = if (inherits(theory, "try-error")) NULL else theory
  )
  class(out) <- "ncar_fit"
  out
}

#' @export
print.ncar_fit <- function(x, ...) {
  cat("Noncausal AR(", x$p, ") with causal GARCH(1,1): three-step fit\n", sep = "")
  cat("Effective sample size: ", length(x$residuals), "\n", sep = "")
  cat("Distribution: ", x$distribution, "\n", sep = "")
  cat("Step 1 phi: ", paste(format(x$step1$phi, digits = 5), collapse = ", "), "\n", sep = "")
  cat("Step 3 phi: ", paste(format(x$step3$phi, digits = 5), collapse = ", "), "\n", sep = "")
  cat("GARCH omega, alpha1, beta1: ",
      paste(format(x$step2$coef[c("omega", "alpha1", "beta1")], digits = 5), collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' @export
summary.ncar_fit <- function(object, ...) {
  garch_coef <- object$step2$coef
  garch_se <- object$step2$se
  out <- list(
    model = paste0("noncausal AR(", object$p, ")-causal GARCH(1,1)"),
    distribution = object$distribution,
    step1 = .ncar_step_table(object$step1$phi, object$step1$se),
    step2 = .ncar_step_table(garch_coef, garch_se),
    step3 = .ncar_step_table(object$step3$phi, object$step3$se),
    convergence = object$step2$convergence,
    settings = object$settings
  )
  class(out) <- "summary.ncar_fit"
  out
}

#' @export
coef.ncar_fit <- function(object, step = c("step3", "step1", "step2", "all"), ...) {
  step <- match.arg(step)
  if (step == "step1") return(object$step1$phi)
  if (step == "step2") return(object$step2$coef)
  if (step == "step3") return(object$step3$phi)
  c(step1 = object$step1$phi, step2 = object$step2$coef, step3 = object$step3$phi)
}

#' @export
vcov.ncar_fit <- function(object, step = c("step3", "step1"), ...) {
  step <- match.arg(step)
  if (step == "step1") return(object$step1$vcov)
  object$step3$vcov
}

#' @export
residuals.ncar_fit <- function(object,
                               type = c("epsilon", "standardized"),
                               ...) {
  type <- match.arg(type)
  if (type == "epsilon") return(object$residuals)
  object$std_residuals
}

#' @export
sigma.ncar_fit <- function(object,
                           type = c("epsilon", "x"),
                           variance = FALSE,
                           ...) {
  type <- match.arg(type)
  if (type == "epsilon") return(if (variance) object$var_epsilon else object$sigma_epsilon)
  if (variance) object$var_x else object$sigma_x
}

#' @export
fitted.ncar_fit <- function(object, step = c("step3", "step1"), ...) {
  step <- match.arg(step)
  if (step == "step1") return(object$step1$fitted)
  object$step3$fitted
}
