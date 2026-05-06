# Compute the modified Qm statistic for m = 1, ..., 20
library(ncargarch)

sim <- ncar_simulate_ar1(
  n = 1500,
  phi = 0.7,
  omega = 0.01,
  alpha = 0.08,
  beta = 0.90,
  dist = "std",
  df = 8,
  seed = 123
)

fit <- ncar_fit(sim$x, p = 1, dist = "norm")
qm <- ncar_test(fit, M = 20, type = "qm")
print(qm)
