# Simulate a noncausal AR(1)-causal GARCH(1,1) process
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
