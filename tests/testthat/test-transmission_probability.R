
# Shared parameters (same values used in test-fsdistr.R)
infectious_period <- 0.9
s0 <- 6
i0 <- 1
beta <- 1.4

N <- s0 + i0

# ip_model = "constant" --------------------------------------------------------

test_that("transmission_probability returns a numeric value (constant)", {
  res <- transmission_probability(beta = beta, n = N, ip_model = "constant",
                                  ip_params = 1 / infectious_period)
  expect_type(res, "double")
})

test_that("transmission_probability is between 0 and 1 (constant)", {
  res <- transmission_probability(beta = beta, n = N, ip_model = "constant",
                                  ip_params = 1 / infectious_period)
  expect_gte(res, 0)
  expect_lte(res, 1)
})

test_that("transmission_probability does not error (constant)", {
  expect_no_error(
    transmission_probability(beta = beta, n = N, ip_model = "constant",
                             ip_params = 1 / infectious_period)
  )
})


# ip_model = "exponential" -----------------------------------------------------

test_that("transmission_probability returns a numeric value (exponential)", {
  res <- transmission_probability(beta = beta, n = N, ip_model = "exponential",
                                  ip_params = 1 / infectious_period)
  expect_type(res, "double")
})

test_that("transmission_probability is between 0 and 1 (exponential)", {
  res <- transmission_probability(beta = beta, n = N, ip_model = "exponential",
                                  ip_params = 1 / infectious_period)
  expect_gte(res, 0)
  expect_lte(res, 1)
})

test_that("transmission_probability does not error (exponential)", {
  expect_no_error(
    transmission_probability(beta = beta, n = N, ip_model = "exponential",
                             ip_params = 1 / infectious_period)
  )
})


# ip_model = "gamma" -----------------------------------------------------------

test_that("transmission_probability returns a numeric value (gamma)", {
  res <- transmission_probability(beta = beta, n = N, ip_model = "gamma",
                                  ip_params = c(2, 2 / infectious_period))
  expect_type(res, "double")
})

test_that("transmission_probability is between 0 and 1 (gamma)", {
  res <- transmission_probability(beta = beta, n = N, ip_model = "gamma",
                                  ip_params = c(2, 2 / infectious_period))
  expect_gte(res, 0)
  expect_lte(res, 1)
})

test_that("transmission_probability does not error (gamma)", {
  expect_no_error(
    transmission_probability(beta = beta, n = N, ip_model = "gamma",
                             ip_params = c(2, 2 / infectious_period))
  )
})


# Vectorised input -------------------------------------------------------------

test_that("transmission_probability handles vectorised beta and s0", {
  res <- transmission_probability(beta = c(beta, 2 * beta), n = c(N, N),
                                  ip_model = "constant",
                                  ip_params = 1 / infectious_period)
  expect_type(res, "double")
  expect_length(res, 2)
  expect_true(all(res >= 0 & res <= 1))
})
