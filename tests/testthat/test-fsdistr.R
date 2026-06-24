


infectious_period <- 0.9
s0 <- 6
i0 <- 1
beta <- 1.4


# fsdistr: ip_model = "constant"---------------------------------------------------------------

test_that("fsdistr returns a vector of length s0 + 1", {
  res <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "constant",
                 ip_params = 1 / infectious_period)
  expect_length(res, s0 + 1)
})

test_that("fsdistr probabilities sum to 1", {
  res <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "constant",
                 ip_params = 1 / infectious_period)
  expect_equal(sum(res), 1, tolerance = 1e-6)
})

test_that("fsdistr probabilities are non-negative", {
  res <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "constant",
                 ip_params = 1 / infectious_period)
  expect_gte(min(res), 0)
})

test_that("fsdistr_internal xmax argument truncates the output", {
  fs0 <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                          ip_params = 1 / infectious_period, ip_model = "constant",
                          xmax = 0)
  fs1 <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                          ip_params = 1 / infectious_period, ip_model = "constant",
                          xmax = 1)
  fsInf <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                            ip_params = 1 / infectious_period, ip_model = "constant",
                            xmax = Inf)
  fs1 <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                          ip_params = 1 / infectious_period, ip_model = "constant",
                          xmax = 1)
  fsInf <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                            ip_params = 1 / infectious_period, ip_model = "constant",
                            xmax = Inf)

  expect_length(fs0, 1)
  expect_length(fs1, 2)
  expect_equal(fsInf[1], fs0[[1]])
  expect_equal(fs1[1], fs0[[1]])
  expect_equal(fsInf[2], fs1[2])
})

test_that("fsdistr equals fsdistr_internal with xmax = Inf", {
  res_full <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "constant",
                      ip_params = 1 / infectious_period)
  res_internal <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                                   ip_params = 1 / infectious_period,
                                   ip_model = "constant", xmax = Inf)
  expect_equal(res_full, res_internal)
})

test_that("fsdistr validates ip_model", {
  expect_snapshot(
    fsdistr(s0 = 5, i0 = 1, beta = 1.5, ip_model = "invalid", ip_params = 1),
    error = TRUE
  )
})

test_that("fsdistr validates ip_params type", {
  expect_snapshot(
    fsdistr(s0 = 5, i0 = 1, beta = 1.5, ip_model = "constant",
            ip_params = "bad"),
    error = TRUE
  )
})


# fsdistr: ip_model = "exponential" -------------------------------------

test_that("fsdistr returns a vector of length s0 + 1 (exponential)", {
  res <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "exponential",
                 ip_params = 1 / infectious_period)
  expect_length(res, s0 + 1)
})

test_that("fsdistr probabilities sum to 1 (exponential)", {
  res <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "exponential",
                 ip_params = 1 / infectious_period)
  expect_equal(sum(res), 1, tolerance = 1e-6)
})

test_that("fsdistr probabilities are non-negative (exponential)", {
  res <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "exponential",
                 ip_params = 1 / infectious_period)
  expect_gte(min(res), 0)
})

test_that("fsdistr_internal xmax argument truncates the output (exponential)", {
  fs0 <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                          ip_params = 1 / infectious_period, ip_model = "exponential",
                          xmax = 0)
  fs1 <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                          ip_params = 1 / infectious_period, ip_model = "exponential",
                          xmax = 1)
  fsInf <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                            ip_params = 1 / infectious_period, ip_model = "exponential",
                            xmax = Inf)

  expect_length(fs0, 1)
  expect_length(fs1, 2)
  expect_equal(fsInf[1], fs0[[1]])
  expect_equal(fs1[1], fs0[[1]])
  expect_equal(fsInf[2], fs1[2])
})

test_that("fsdistr equals fsdistr_internal with xmax = Inf (exponential)", {
  res_full <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "exponential",
                      ip_params = 1 / infectious_period)
  res_internal <- fsdistr_internal(s0 = s0, i0 = i0, beta = beta,
                                   ip_params = 1 / infectious_period,
                                   ip_model = "exponential", xmax = Inf)
  expect_equal(res_full, res_internal)
})


# fsdistr: ip_model = "gamma" -------------------------------------------

test_that("fsdistr returns a vector of length s0 + 1 (gamma)", {
  res <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "gamma",
                 ip_params = c(2, 2 / infectious_period))
  expect_length(res, s0 + 1)
})

test_that("fsdistr probabilities sum to 1 (gamma)", {
  res <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "gamma",
                 ip_params = c(2, 2 / infectious_period))
  expect_equal(sum(res), 1, tolerance = 1e-6)
})

test_that("fsdistr probabilities are non-negative (gamma)", {
  res <- fsdistr(s0 = s0, i0 = i0, beta = beta, ip_model = "gamma",
                 ip_params = c(2, 2 / infectious_period))
  expect_gte(min(res), 0)
})


# fsdistr_mt ------------------------------------------------------------

test_that("fsdistr_mt single-group matches fsdistr (constant)", {

  res_mt <- fsdistr_mt(s0 = s0, i0 = i0, beta = matrix(beta),
                       ip_model = "constant",
                       ip_params = list(1 / infectious_period))

  res_1 <- fsdistr(s0 = s0, i0 = i0, beta = beta,
                   ip_model = "constant",
                   ip_params = 1 / infectious_period)

  expect_equal(unname(res_mt), res_1, tolerance = 1e-6)

})

test_that("fsdistr_mt single-group matches fsdistr (expoential)", {

  res_mt <- fsdistr_mt(s0 = s0, i0 = i0, beta = matrix(beta),
                       ip_model = "exponential",
                       ip_params = list(1 / infectious_period))

  res_1 <- fsdistr(s0 = s0, i0 = i0, beta = beta,
                   ip_model = "exponential",
                   ip_params = 1 / infectious_period)

  expect_equal(unname(res_mt), res_1, tolerance = 1e-6)
})

test_that("fsdistr_mt single-group matches fsdistr (gamma)", {

  res_mt <- fsdistr_mt(s0 = s0, i0 = i0, beta = matrix(beta),
                       ip_model = "gamma",
                       ip_params = list(c(1, 1 / infectious_period)))

  res_1 <- fsdistr(s0 = s0, i0 = i0, beta = beta,
                   ip_model = "gamma",
                   ip_params = c(1, 1 / infectious_period))

  expect_equal(unname(res_mt), res_1, tolerance = 1e-6)
})

test_that("fsdistr_mt probabilities sum to 1", {
  beta2 <- matrix(c(1.5, 0.5, 0.5, 1.5), nrow = 2)
  res <- fsdistr_mt(s0 = c(4, 4), i0 = c(1, 0), beta = beta2)
  expect_equal(sum(res), 1, tolerance = 1e-4)
})

test_that("fsdistr_mt probabilities are non-negative", {
  beta2 <- matrix(c(1.5, 0.5, 0.5, 1.5), nrow = 2)
  res <- fsdistr_mt(s0 = c(4, 4), i0 = c(1, 0), beta = beta2)
  expect_gte(min(res), 0)
})

test_that("fsdistr_mt result names encode outcomes correctly", {
  res <- fsdistr_mt(s0 = c(2, 2), i0 = c(1, 0), beta = matrix(c(1, 0.5, 0.5, 1), nrow = 2))
  expect_match(names(res)[1], "^\\d+-\\d+$")
})


test_that("fsdistr_mt result marignal distribution matches fsdistr", {

  s0 <- c(3, 1)
  i0 <- c(1,1)

  s0_tot <- sum(s0)
  i0_tot <- sum(i0)


  beta_coef <- 1.5
  beta_mat <- beta_coef * matrix(c(1.0, 1.0, 1.0, 1.0), nrow = 2, ncol = 2, byrow = TRUE)


  fsdistr_res <- fsdistr(s0 = s0_tot, i0 = i0_tot, beta = 1.5, ip_model = 'constant')

  # Compute mt distribution, then marginalize to get the overall distribution.
  fsdistr_mt_res <- fsdistr_mt(i0 = i0, s0 = s0, beta = beta_mat, ip_model = 'constant', return_df = TRUE)
  total_i <- rowSums(fsdistr_mt_res[,-ncol(fsdistr_mt_res), drop=FALSE])
  fsdistr_mt_res_tot <- sapply(split(x = fsdistr_mt_res$probability, f = total_i), FUN = sum)

  # also compute the mt distributin with reverse ordering of the groups
  fsdistr_mt_res_rev <- fsdistr_mt(i0 = rev(i0), s0 = rev(s0), beta = beta_mat, ip_model = 'constant', return_df = TRUE)
  total_i_rev <- rowSums(fsdistr_mt_res_rev[,-ncol(fsdistr_mt_res_rev), drop=FALSE])
  fsdistr_mt_res_tot_rev <- sapply(split(x = fsdistr_mt_res_rev$probability, f = total_i_rev), FUN = sum)



  expect_equal(unname(fsdistr_mt_res_tot), fsdistr_res, tolerance = 1e-6)
  expect_equal(unname(fsdistr_mt_res_tot_rev), fsdistr_res, tolerance = 1e-6)

})




# make_multitype_state_table ----
# It is important that the ordering of the outcomes made by make_multitype_state_table() is
# correct and does not change in the future (it relies on the built-in expand.grid function).

test_that("make_multitype_state_table gives desired result", {

  sm <- make_multitype_state_table(c(3,3))

  # Hard coded desired output
  sm_c1 <- c(0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3)
  sm_c2 <- c(0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)

  expect_equal(sm[,1], sm_c1)
  expect_equal(sm[,2], sm_c2)

})






