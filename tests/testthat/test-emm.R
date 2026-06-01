

s0 <- 20
i0 <- 2
beta <- 1.6
gamma <- 1.1
beta0 <- 0.2


# make_markov_state_table -------------------------------------------------------

test_that("make_markov_state_table returns correct number of states", {
  N <- s0 + i0
  st <- make_markov_state_table(N)
  expect_equal(nrow(st), sum(1:(N + 1)))
})

test_that("make_markov_state_table columns are S and I counts", {
  N <- 3
  st <- make_markov_state_table(N)
  expect_equal(ncol(st), 2)
  # Every row must satisfy S >= 0, I >= 0, S + I <= N
  expect_equal(min(st[, 1]), 0)
  expect_equal(min(st[, 2]), 0)
  expect_equal(max(rowSums(st)), N)
})


# make_transition_matrix -------------------------------------------------

test_that("make_transition_matrix dense rows sum to 1 (probabilities)", {
  st <- make_markov_state_table(s0 + i0)
  tm <- make_transition_matrix(st, s0 = s0, beta = beta, gamma = gamma, sparse = FALSE)
  expect_equal(rowSums(tm), rep(1, nrow(tm)), tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("make_transition_matrix sparse rows sum to 1 (probabilities)", {
  st <- make_markov_state_table(s0 + i0)
  tm <- make_transition_matrix(st, s0 = s0, beta = beta, gamma = gamma, sparse = TRUE)
  expect_equal(as.numeric(Matrix::rowSums(tm)), rep(1, nrow(tm)), tolerance = 1e-6)
})

test_that("make_transition_matrix sparse and dense are equivalent (probabilities)", {
  st <- make_markov_state_table(s0 + i0)
  tm_dense <- make_transition_matrix(st, s0 = s0, beta = beta, gamma = gamma, sparse = FALSE)
  tm_sparse <- make_transition_matrix(st, s0 = s0, beta = beta, gamma = gamma, sparse = TRUE)

  expect_equal(dim(tm_dense), dim(tm_sparse))
  expect_equal(colnames(tm_dense), colnames(tm_sparse))
  expect_equal(rownames(tm_dense), rownames(tm_sparse))
  expect_equal(as.matrix(tm_sparse), tm_dense, tolerance = 1e-10)
})


test_that("make_transition_matrix sparse and dense are equivalent (rates)", {
  st <- make_markov_state_table(s0 + i0)
  tm_dense <- make_transition_matrix(st, s0 = s0, beta = beta, gamma = gamma, sparse = FALSE, elements = 'rates')
  tm_sparse <- make_transition_matrix(st, s0 = s0, beta = beta, gamma = gamma, sparse = TRUE, elements = 'rates')

  expect_equal(dim(tm_dense), dim(tm_sparse))
  expect_equal(colnames(tm_dense), colnames(tm_sparse))
  expect_equal(rownames(tm_dense), rownames(tm_sparse))
  expect_equal(as.matrix(tm_sparse), tm_dense, tolerance = 1e-10)
})


# emmdt ------------------------------------------------------------------

test_that("emmdt fs_distr has length s0 + 1", {
  res <- emmdt(s0 = s0, i0 = i0, beta = beta, gamma = gamma)
  expect_length(res$fs_distr, s0 + 1)
})

test_that("emmdt fs_distr sums to 1", {
  res <- emmdt(s0 = s0, i0 = i0, beta = beta, gamma = gamma)
  expect_equal(sum(res$fs_distr), 1, tolerance = 1e-6)
})

test_that("emmdt fs_distr probabilities are non-negative", {
  res <- emmdt(s0 = s0, i0 = i0, beta = beta, gamma = gamma)
  expect_gte(min(res$fs_distr), 0)
})

test_that("emmdt sparse and dense give same fs_distr", {
  res_sparse <- emmdt(s0 = s0, i0 = i0, beta = beta, gamma = gamma, sparse = TRUE)
  res_dense <- emmdt(s0 = s0, i0 = i0, beta = beta, gamma = gamma, sparse = FALSE)
  expect_equal(res_sparse$fs_distr, res_dense$fs_distr, tolerance = 1e-6)
})

test_that("emmdt returns expected list elements", {
  res <- emmdt(s0 = s0, i0 = i0, beta = beta, gamma = gamma)
  expect_named(res, c("transition_matrix", "qmat", "rmat", "fmat", "solution_mat", "fs_distr"))
})

test_that("emmdt solution_mat rows sum to 1", {
  res <- emmdt(s0 = s0, i0 = i0, beta = beta, gamma = gamma)
  expect_equal(
    as.numeric(Matrix::rowSums(res$solution_mat)),
    rep(1, nrow(res$solution_mat)),
    tolerance = 1e-6
  )
})

test_that("emmdt requires i0 >= 1", {
  expect_snapshot(emmdt(s0 = 10, i0 = 0, beta = 1.5, gamma = 1.0), error = TRUE)
})

test_that("emmdt requires s0 >= 1", {
  expect_snapshot(emmdt(s0 = 0, i0 = 1, beta = 1.5, gamma = 1.0), error = TRUE)
})


test_that("emmdt matches  fsdist (exponential)", {
  res_emmdt <- emmdt(s0 = s0, i0 = i0, beta = beta, gamma = gamma)

  res_fsdistr <- fsdistr(s0 = s0, i0 = i0, beta = beta,
                         ip_model = "exponential",
                         ip_params = gamma)

  expect_true(all(abs(res_emmdt$fs_distr - res_fsdistr) <= 1e-10))
})



# emmct ------------------------------------------------------------------

test_that("emmct fs_distr has length s0 + 1", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_length(res1$fs_distr, s0 + 1)
  expect_length(res2$fs_distr, s0 + 1)
  expect_length(res3$fs_distr, s0 + 1)
})

test_that("emmct fs_distr sums to 1", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_equal(sum(res1$fs_distr), 1, tolerance = 1e-6)
  expect_equal(sum(res2$fs_distr), 1, tolerance = 1e-6)
  expect_equal(sum(res3$fs_distr), 1, tolerance = 1e-6)

})

test_that("emmct fs_distr probabilities are non-negative", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_gte(min(res1$fs_distr), 0)
  expect_gte(min(res2$fs_distr), 0)
  expect_gte(min(res3$fs_distr), 0)

})

test_that("emmct transition matrix rows sum to 1", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_all_true(abs(Matrix::rowSums(res1$transition_matrix)-1) <= 1e-7)
  expect_all_true(abs(Matrix::rowSums(res2$transition_matrix)-1) <= 1e-7)
  expect_all_true(abs(Matrix::rowSums(res2$transition_matrix)-1) <= 1e-7)
})

test_that("emmct transition matrix probabilities are non-negative", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_gte(min(res1$transition_matrix), 0)
  expect_gte(min(res2$transition_matrix), 0)
  expect_gte(min(res3$transition_matrix), 0)
})


test_that("emmct sparse and dense are equivalent", {
  res_s <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.1, sparse = TRUE)
  res_d <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.1, sparse = FALSE)

  expect_equal(res_s$fs_distr, res_d$fs_distr)
  expect_equal(as.matrix(res_s$transition_matrix), res_d$transition_matrix)

})


test_that("emmct fs_distr has length s0 + 1 (beta0 != 0)", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_length(res1$fs_distr, s0 + 1)
  expect_length(res2$fs_distr, s0 + 1)
  expect_length(res3$fs_distr, s0 + 1)
})

test_that("emmct fs_distr sums to 1 (beta0 != 0)", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_equal(sum(res1$fs_distr), 1, tolerance = 1e-6)
  expect_equal(sum(res2$fs_distr), 1, tolerance = 1e-6)
  expect_equal(sum(res3$fs_distr), 1, tolerance = 1e-6)

})

test_that("emmct fs_distr probabilities are non-negative (beta0 != 0)", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_gte(min(res1$fs_distr), 0)
  expect_gte(min(res2$fs_distr), 0)
  expect_gte(min(res3$fs_distr), 0)

})

test_that("emmct transition matrix rows sum to 1 (beta0 != 0)", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_all_true(abs(Matrix::rowSums(res1$transition_matrix)-1) <= 1e-7)
  expect_all_true(abs(Matrix::rowSums(res2$transition_matrix)-1) <= 1e-7)
  expect_all_true(abs(Matrix::rowSums(res2$transition_matrix)-1) <= 1e-7)
})

test_that("emmct transition matrix probabilities are non-negative (beta0 != 0)", {
  res1 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_gte(min(res1$transition_matrix), 0)
  expect_gte(min(res2$transition_matrix), 0)
  expect_gte(min(res3$transition_matrix), 0)
})


test_that("emmct sparse and dense are equivalent (beta0 != 0)", {
  res_s <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.1, sparse = TRUE)
  res_d <- emmct(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.1, sparse = FALSE)

  expect_equal(res_s$fs_distr, res_d$fs_distr)
  expect_equal(as.matrix(res_s$transition_matrix), res_d$transition_matrix)

})

# emmct_sis ------------------------------------------------------------------

test_that("emmct_sis sparse and dense are equivalent", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.1, sparse = TRUE)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.1, sparse = FALSE)

  expect_equal(res1$state_prob,res2$state_prob)
})

test_that("emmct_sis state_prob has length s0 + 1", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_length(res1$state_prob, s0 + i0 + 1)
  expect_length(res2$state_prob, s0 + i0 + 1)
  expect_length(res3$state_prob, s0 + i0 + 1)
})

test_that("emmct_sis state_prob sums to 1", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_equal(sum(res1$state_prob), 1, tolerance = 1e-6)
  expect_equal(sum(res2$state_prob), 1, tolerance = 1e-6)
  expect_equal(sum(res3$state_prob), 1, tolerance = 1e-6)

})

test_that("emmct_sis state_prob probabilities are non-negative", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_gte(min(res1$state_prob), 0)
  expect_gte(min(res2$state_prob), 0)
  expect_gte(min(res3$state_prob), 0)

})

test_that("emmct_sis transition matrix rows sum to 1", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_all_true(abs(Matrix::rowSums(res1$transition_matrix)-1) <= 1e-7)
  expect_all_true(abs(Matrix::rowSums(res2$transition_matrix)-1) <= 1e-7)
  expect_all_true(abs(Matrix::rowSums(res2$transition_matrix)-1) <= 1e-7)
})

test_that("emmct_sis transition matrix probabilities are non-negative", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, time = 5.5)

  expect_gte(min(res1$transition_matrix), 0)
  expect_gte(min(res2$transition_matrix), 0)
  expect_gte(min(res3$transition_matrix), 0)
})

### ### ### ### ### ### ### ### ###

test_that("emmct_sis sparse and dense are equivalent (beta0 != 0)", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.1, sparse = TRUE)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.1, sparse = FALSE)

  expect_equal(res1$state_prob,res2$state_prob)
})

test_that("emmct_sis state_prob has length s0 + 1 (beta0 != 0)", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_length(res1$state_prob, s0 + i0 + 1)
  expect_length(res2$state_prob, s0 + i0 + 1)
  expect_length(res3$state_prob, s0 + i0 + 1)
})

test_that("emmct_sis state_prob sums to 1 (beta0 != 0)", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_equal(sum(res1$state_prob), 1, tolerance = 1e-6)
  expect_equal(sum(res2$state_prob), 1, tolerance = 1e-6)
  expect_equal(sum(res3$state_prob), 1, tolerance = 1e-6)

})

test_that("emmct_sis state_prob probabilities are non-negative (beta0 != 0)", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_gte(min(res1$state_prob), 0)
  expect_gte(min(res2$state_prob), 0)
  expect_gte(min(res3$state_prob), 0)

})

test_that("emmct_sis transition matrix rows sum to 1 (beta0 != 0)", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_all_true(abs(Matrix::rowSums(res1$transition_matrix)-1) <= 1e-7)
  expect_all_true(abs(Matrix::rowSums(res2$transition_matrix)-1) <= 1e-7)
  expect_all_true(abs(Matrix::rowSums(res2$transition_matrix)-1) <= 1e-7)
})

test_that("emmct_sis transition matrix probabilities are non-negative (beta0 != 0)", {
  res1 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 0.1)
  res2 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 1.2)
  res3 <- emmct_sis(s0 = s0, i0 = i0, beta = beta, gamma = gamma, beta0 = beta0, time = 5.5)

  expect_gte(min(res1$transition_matrix), 0)
  expect_gte(min(res2$transition_matrix), 0)
  expect_gte(min(res3$transition_matrix), 0)
})










