.ncar_extract_column <- function(data, column = "Adjusted") {
  column <- as.character(column)[1L]

  if (is.numeric(data) && is.null(dim(data))) {
    return(list(values = as.numeric(data), dates = NULL, raw = data, column = column))
  }

  dates <- NULL
  if (requireNamespace("zoo", quietly = TRUE) && zoo::is.zoo(data)) {
    dates <- zoo::index(data)
  }

  cn <- colnames(data)
  if (is.null(cn)) {
    values <- as.numeric(data)
    return(list(values = values, dates = dates, raw = data, column = column))
  }

  target <- switch(tolower(column),
    adjusted = "adjusted",
    adj = "adjusted",
    close = "close",
    cl = "close",
    open = "open",
    high = "high",
    low = "low",
    volume = "volume",
    tolower(column)
  )

  idx <- grep(paste0("(^|\\.)", target, "$"), tolower(cn))
  if (length(idx) == 0L) {
    idx <- which(tolower(cn) == target)
  }
  if (length(idx) == 0L) {
    stop("Could not find column '", column, "'. Available columns: ", paste(cn, collapse = ", "), call. = FALSE)
  }
  idx <- idx[1L]

  list(values = as.numeric(data[, idx]), dates = dates, raw = data, column = cn[idx])
}

#' Download market data using quantmod
#'
#' Convenience wrapper around `quantmod::getSymbols()` with
#' `auto.assign = FALSE`.
#'
#' @param symbol Symbol to download, for example `"XRP-USD"`.
#' @param from Start date.
#' @param to End date.
#' @param src Data source passed to `quantmod::getSymbols()`.
#' @param ... Additional arguments passed to `quantmod::getSymbols()`.
#'
#' @return Usually an xts object, depending on the source.
#' @export
ncar_get_data <- function(symbol,
                          from,
                          to = Sys.Date(),
                          src = "yahoo",
                          ...) {
  .ncar_require_namespace("quantmod")
  quantmod::getSymbols(
    Symbols = symbol,
    src = src,
    from = from,
    to = to,
    auto.assign = FALSE,
    ...
  )
}

#' Prepare downloaded or supplied data for `ncar_fit()`
#'
#' @param data Numeric vector, xts/zoo object, or data frame.
#' @param column Column to use when `data` has columns. Common choices are
#'   `"Adjusted"` and `"Close"`.
#' @param transform Transformation: `"price"`, `"log_price"`, `"return"`, or
#'   `"detrended_price"`.
#' @param demean Logical. If `TRUE`, subtract the sample mean after the selected
#'   transformation.
#' @param detrend_args List of arguments passed to [ncar_detrend()] when
#'   `transform = "detrended_price"`.
#'
#' @return A list with prepared vector `x`, `dates`, and metadata.
#' @export
ncar_prepare_data <- function(data,
                              column = "Adjusted",
                              transform = c("price", "log_price", "return", "detrended_price"),
                              demean = FALSE,
                              detrend_args = list()) {
  transform <- match.arg(transform)
  extracted <- .ncar_extract_column(data, column = column)
  price <- extracted$values
  dates <- extracted$dates

  ok <- is.finite(price)
  if (!all(ok)) {
    price <- price[ok]
    if (!is.null(dates)) dates <- dates[ok]
  }

  detrend <- NULL
  if (transform == "price") {
    x <- price
  } else if (transform == "log_price") {
    if (any(price <= 0)) stop("log_price requires strictly positive prices.", call. = FALSE)
    x <- log(price)
  } else if (transform == "return") {
    if (any(price <= 0)) stop("return requires strictly positive prices.", call. = FALSE)
    x <- diff(log(price))
    if (!is.null(dates)) dates <- dates[-1L]
  } else {
    detrend <- do.call(ncar_detrend, c(list(x = price), detrend_args))
    x <- detrend$detrended
  }

  x_mean <- if (demean) mean(x) else 0
  x <- x - x_mean

  out <- list(
    x = as.numeric(x),
    dates = dates,
    raw = extracted$raw,
    price = price,
    column = extracted$column,
    transform = transform,
    demean = demean,
    mean_removed = x_mean,
    detrend = detrend
  )
  class(out) <- "ncar_prepared_data"
  out
}

#' Download, prepare, and fit one symbol
#'
#' @param symbol Symbol to download.
#' @param from Start date.
#' @param to End date.
#' @param src Data source passed to [ncar_get_data()].
#' @param column Price column passed to [ncar_prepare_data()].
#' @param transform Transformation passed to [ncar_prepare_data()].
#' @param p Noncausal AR order passed to [ncar_fit()].
#' @param dist Conditional distribution passed to [ncar_fit()].
#' @param demean Logical passed to [ncar_prepare_data()].
#' @param detrend_args List passed to [ncar_prepare_data()].
#' @param ... Additional arguments passed to [ncar_fit()].
#'
#' @return An `ncar_fit` object with additional `download` and `prepared_data`
#'   components.
#' @export
ncar_fit_symbol <- function(symbol,
                            from,
                            to = Sys.Date(),
                            src = "yahoo",
                            column = "Adjusted",
                            transform = c("detrended_price", "price", "log_price", "return"),
                            p = 1L,
                            dist = c("norm", "std", "sstd"),
                            demean = FALSE,
                            detrend_args = list(),
                            ...) {
  transform <- match.arg(transform)
  dist <- .ncar_match_dist(dist)
  raw <- ncar_get_data(symbol = symbol, from = from, to = to, src = src)
  prepared <- ncar_prepare_data(
    data = raw,
    column = column,
    transform = transform,
    demean = demean,
    detrend_args = detrend_args
  )
  fit <- ncar_fit(
    x = prepared$x,
    p = p,
    dist = dist,
    demean = FALSE,
    dates = prepared$dates,
    ...
  )
  fit$download <- list(symbol = symbol, from = from, to = to, src = src, column = column)
  fit$prepared_data <- prepared
  fit
}
