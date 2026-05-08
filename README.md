# ncargarch
Noncausal AR models with causal GARCH volatility

`ncargarch` contains R functions for simulation, estimation, diagnostic testing, detrending, and empirical workflows for noncausal autoregressive models driven by causal GARCH(1,1) volatility.

The core model is

```text
X_t = phi_1 X_{t+1} + ... + phi_p X_{t+p} + epsilon_t

epsilon_t = sigma_t eta_t
sigma_t^2 = omega + alpha epsilon_{t-1}^2 + beta sigma_{t-1}^2
```
## Related working paper

This repository provides R code associated with the working paper:

> Velásquez-Gaviria, D. and Zakoïan, J.-M. (2026).  
> **Noncausal AR processes driven by causal GARCH volatility**.  
> CREST Working Papers Series No. 2026-02, Center for Research in Economics and Statistics.  
> RePEc handle: `RePEc:crs:wpaper:2026-02`.

Working paper link: https://crest.science/wp-content/uploads/2026/02/2026-02.pdf


## Main functions

| Function | Purpose |
|---|---|
| `ncar_simulate_ar1()` | Simulate AR(1)-causal GARCH(1,1) |
| `ncar_simulate_arp()` | Simulate AR(p)-causal GARCH(1,1) |
| `ncar_fit()` | Three-step LS/QML/WLS estimation |
| `ncar_test()` | Portmanteau tests; default `type = "qm"` |
| `ncar_components()` | Extract residuals, standardized residuals, volatility, and fitted values |
| `ncar_export()` | Export estimated components to CSV or RDS |
| `ncar_detrend()` | One-sided local-level detrending with ADF selection |
| `ncar_get_data()` | Download data using `quantmod` |
| `ncar_prepare_data()` | Clean and transform downloaded data |
| `ncar_fit_symbol()` | Download, prepare, and fit one symbol |

## Installation from GitHub

GitHub:

```r
# install.packages("remotes")
remotes::install_github("dvelas18/ncargarch")
```

For local development:

```r
# install.packages("devtools")
devtools::load_all()
```

Suggested packages:

```r
install.packages(c("rugarch", "quantmod", "xts", "zoo", "urca", "testthat"))
```

## Simulate AR(1)-GARCH(1,1)

```r
library(ncargarch)

sim <- ncar_simulate_ar1(
  n = 1000,
  phi = 0.9,
  omega = 0.01,
  alpha = 0.10,
  beta = 0.85,
  dist = "std",
  df = 10,
  seed = 123
)

print(sim)
plot(sim)
```

## Simulate with an imposed shock

```r
sim_shock <- ncar_simulate_ar1(
  n = 700,
  phi = 0.9,
  omega = 0.1,
  alpha = 0.10,
  beta = 0.75,
  dist = "std",
  df = 10,
  shock = list(time = 500, value = 40, type = "eta"),
  seed = 123
)

plot(sim_shock, which = c("x", "epsilon_volatility", "x_volatility"))
```

## Three-step estimation

```r
fit <- ncar_fit(sim$x, p = 1, dist = "norm", se = "theory")
summary(fit)
```

Extract components:

```r
comp <- ncar_components(fit)
head(comp)
```

Available components include:

- `epsilon_hat`: estimated noncausal AR residuals;
- `eta_hat`: standardized residuals;
- `sigma_epsilon_hat`: GARCH innovation volatility;
- `var_epsilon_hat`: GARCH innovation variance;
- `sigma_x_hat`: AR(1) observed-process volatility under the innovation filtration, when available;
- `var_x_hat`: corresponding variance measure.

## Qm diagnostic test

```r
qm <- ncar_test(fit, M = 20, type = "qm")
print(qm)
```

Benchmark diagnostics are also available:

```r
ncar_test(fit, m = c(5, 10, 20), type = "ljungbox")
ncar_test(fit, m = c(5, 10, 20), type = "boxpierce")
```

The exact modified `type = "qm"` implementation is currently restricted to AR(1), matching the current theoretical result.

## Detrending workflow

```r
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

fit_sol <- ncar_fit(prep$x, p = 1, dist = "norm", dates = prep$dates, se = "theory")
summary(fit_sol)
ncar_test(fit_sol, m = c(1, 5, 10, 15, 20), type = "qm")
ncar_plot_volatility(fit_sol)
```

