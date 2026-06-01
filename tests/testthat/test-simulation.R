# Simulation-based tests: compare fsdistr() / fsdistr_mt() against Monte Carlo results.
#
# These tests are stochastic but use fixed seeds for reproducibility.
# They are skipped on CRAN because they are slower than typical unit tests.
#
# Tolerances are chosen to give a comfortable margin above the Monte Carlo
# standard error at the chosen nsim (each pointwise SE ≈ sqrt(p(1-p)/nsim)).


# fsdistr() – single-type -----------------------------------------------------

test_that("fsdistr matches simulation (exponential infectious period)", {
  skip_on_cran()

  i0   <- 2
  s0   <- 5
  beta <- 1.5
  nsim <- 15000

  set.seed(4812)
  sim <- simulation_distr(
    nsim      = nsim,
    i0        = i0,
    s0        = matrix(s0),
    beta      = matrix(beta),
    ip_model  = "exponential",
    ip_params = list(0.91)
  )

  analytic <- fsdistr(s0 = s0, i0 = i0, beta = beta,
                      ip_model = "exponential", ip_params = 0.91)

  expect_true(
    all(abs(sim - analytic) <= 0.02),
    label = paste0(
      "max deviation = ", round(max(abs(sim - analytic)), 4),
      " (tolerance 0.02)"
    )
  )
})


test_that("fsdistr matches simulation (constant infectious period)", {
  skip_on_cran()

  i0   <- 2
  s0   <- 5
  beta <- 1.5
  nsim <- 15000

  set.seed(7341)
  sim <- simulation_distr(
    nsim      = nsim,
    i0        = i0,
    s0        = matrix(s0),
    beta      = matrix(beta),
    ip_model  = "constant",
    ip_params = list(1.1)
  )

  analytic <- fsdistr(s0 = s0, i0 = i0, beta = beta,
                      ip_model = "constant", ip_params = 1.1)

  expect_true(
    all(abs(sim - analytic) <= 0.02),
    label = paste0(
      "max deviation = ", round(max(abs(sim - analytic)), 4),
      " (tolerance 0.02)"
    )
  )
})


test_that("fsdistr matches simulation (gamma infectious period)", {
  skip_on_cran()

  i0   <- 2
  s0   <- 5
  beta <- 1.5
  nsim <- 15000

  set.seed(2903)
  sim <- simulation_distr(
    nsim      = nsim,
    i0        = i0,
    s0        = matrix(s0),
    beta      = matrix(beta),
    ip_model  = "gamma",
    ip_params = list(c(1.1, 0.9))
  )

  analytic <- fsdistr(s0 = s0, i0 = i0, beta = beta,
                      ip_model = "gamma", ip_params = c(1.1, 0.9))

  expect_true(
    all(abs(sim - analytic) <= 0.02),
    label = paste0(
      "max deviation = ", round(max(abs(sim - analytic)), 4),
      " (tolerance 0.02)"
    )
  )
})


# fsdistr_mt() – multi-type ---------------------------------------------------

test_that("fsdistr_mt matches simulation (3-group, mixed infectious periods)", {
  skip_on_cran()

  s0 <- c(2, 1, 1)
  i0 <- c(1, 0, 0)
  m  <- length(s0)

  my_beta <- matrix(
    c(0.4, 1.2, 0.4,
      1.0, 0.1, 0.1,
      0.5, 0.2, 1.09),
    ncol = 3, nrow = 3, byrow = TRUE
  )

  my_ip_params <- list(c(1.1, 0.9), 0.92, 1.15)
  my_ip_model  <- c("gamma", "constant", "exponential")

  nsim <- 1000000

  set.seed(6157)
  sim_mat <- simulation_distr_mt(
    nsim      = nsim,
    i0        = i0,
    s0        = s0,
    beta      = my_beta,
    ip_model  = my_ip_model,
    ip_params = my_ip_params
  )

  sim_probs <- sim_mat[, m + 1] / nsim

  analytic <- fsdistr_mt(
    s0        = s0,
    i0        = i0,
    beta      = my_beta,
    ip_model  = my_ip_model,
    ip_params = my_ip_params
  )

  expect_true(
    all(abs(sim_probs - analytic) <= 0.03),
    label = paste0(
      "max deviation = ", round(max(abs(sim_probs - analytic)), 4),
      " (tolerance 0.03)"
    )
  )
})
