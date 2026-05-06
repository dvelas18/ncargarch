# Simulate a shock experiment similar to the paper
library(ncargarch)

sim <- ncar_simulate_ar1(
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

plot(sim, which = c("x", "epsilon_volatility", "x_volatility"))
