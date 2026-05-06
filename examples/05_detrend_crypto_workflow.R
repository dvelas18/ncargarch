# Download, detrend, and fit a cryptocurrency price series
# Requires internet access and the suggested packages quantmod, urca, and rugarch.
library(ncargarch)

raw <- ncar_get_data("SOL-USD", from = "2020-04-10", to = "2024-07-01")
prep <- ncar_prepare_data(
  raw,
  column = "Adjusted",
  transform = "detrended_price",
  detrend_args = list(
    r_grid = exp(seq(log(1e-6), log(5), length.out = 120)),
    alpha = 0.01,
    deterministic = "drift"
  )
)

fit <- ncar_fit(prep$x, p = 1, dist = "norm", dates = prep$dates, se = "theory")
summary(fit)
ncar_test(fit, m = c(1, 5, 10, 15, 20), type = "qm")
ncar_plot_volatility(fit)
