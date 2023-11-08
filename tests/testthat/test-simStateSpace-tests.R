## ---- test-simStateSpace-tests
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      paste(text, "positive definite"),
      {
        testthat::expect_true(
          simStateSpace:::.TestPositiveDefinite(
            matrix(
              data = diag(2),
              nrow = 2
            )
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "nonpositive definite"),
      {
        testthat::expect_false(
          simStateSpace:::.TestPositiveDefinite(
            matrix(
              data = c(-2, 1, 1, 0),
              nrow = 2
            )
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "stationary"),
      {
        testthat::expect_true(
          simStateSpace:::.TestStationarity(
            matrix(
              data = c(0.5, 0.3, 0.2, 0.4),
              nrow = 2
            )
          )
        )
      }
    )
    testthat::test_that(
      paste(text, "nonstationary"),
      {
        testthat::expect_false(
          simStateSpace:::.TestStationarity(
            matrix(
              data = c(0.9, -0.5, 0.8, 0.7),
              nrow = 2
            )
          )
        )
      }
    )
  },
  text = "test-simStateSpace-tests"
)
