library(testthat)
library(RAEN)

#test_check("RAEN")
test_that("str_length is number of characters", {
  expect_equal(nchar("a"), 1)
  expect_equal(nchar("ab"), 2)
  expect_equal(nchar("abc"), 3)
})