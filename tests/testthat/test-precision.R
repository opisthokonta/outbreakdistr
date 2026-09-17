



test_that("as_mpfr_matrix output as expected", {

  my_matrix <- 1.001154 * matrix(c(1.2, 0.9, 3.0, 0.2, 9, 4), nrow = 3, ncol = 2, byrow = TRUE)
  my_matrix_precise <- as_mpfr_matrix(my_matrix, prec = 64)

  # The numeric elements should be in the same place in the input and output.
  expect_true(all(abs(my_matrix_precise - my_matrix) <= 0.0000001))

})



test_that("fsdistr works with prec argument", {


  res1 <- fsdistr(s0 = 10, i0 = 1, beta = 1.012)
  res2 <- fsdistr(s0 = 10, i0 = 1, beta = 1.012, prec = 64)

  expect_equal(res1, res2, tolerance = 0.000001)

})








