test_that("errors work", {
  endogenous = list(a = 1, b = 1)
  exogenous = list(x = 2, y = 3, z = 3)

  expect_equal(
    errors(
      formulas = list(a = a ~ x + y, b = b ~ y + z),
      exogenous = exogenous,
      endogenous = endogenous
    ),
    list(a = -4, b = -5)
  )
})
