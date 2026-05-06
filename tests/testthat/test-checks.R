test_that("parameter checker accepts a standard AR(1)-GARCH design", {
  res <- ncar_check_params(0.9, omega = 0.01, alpha = 0.1, beta = 0.8, stop_on_error = FALSE)
  expect_true(res$ok)
  expect_true(all(res$root_modulus > 1))
})

test_that("AR root checker rejects invalid noncausal AR(1)", {
  res <- ncar_check_params(1.1, omega = 0.01, alpha = 0.1, beta = 0.8, stop_on_error = FALSE)
  expect_false(res$ok)
})
