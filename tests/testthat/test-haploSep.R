context("Separates haplotypes")

test_that("response to weird inputs", {
  expect_error(haploSep(0))
  expect_error(haploSep(matrix(NA,5,2)))
  expect_error(haploSep(matrix(Inf,5,2)))
  expect_error(haploSep(matrix(0,5,2),3))
})

test_that("validity of the reconstruction", {
  M    = 3
  n    = 5
  data = matrix(abs(rnorm(M*n)),n,M)
  re   = haploSep(data)
  expect_true(setequal(c(re$haploStr),c(0,1))) # check alphabet structure
  expect_true(all(re$omega >= 0)) # check nonnegativity of weights
  # match of dimension
  expect_equal(attr(re,"nHaplo"), ncol(re$haploStr))
  expect_equal(ncol(re$haploStr), nrow(re$haploFrq))
  expect_equal(nrow(re$haploStr), nrow(data))
  expect_equal(ncol(data), ncol(re$haploFrq))
})
