# Estimate the model using simulated data
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

fit <- ncar_fit(sim$x, p = 1, dist = "norm", se = "theory")
summary(fit)

components <- ncar_components(fit)
head(components)
plot(fit)
