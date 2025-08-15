## ---- test-simStateSpace-spectral-radius
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      paste0(text, "x < 1"),
      {
        testthat::skip_on_cran()
        x <- matrix(
          data = c(
            0.5, 0.3,
            0.2, 0.4
          ),
          nrow = 2
        )
        testthat::expect_true(
          SpectralRadius(x = x) < 1
        )
      }
    )
    testthat::test_that(
      paste0(text, "x > 1"),
      {
        testthat::skip_on_cran()
        x <- matrix(
          data = c(
            1.2, 0.3,
            0.4, 0.9
          ),
          nrow = 2
        )
        testthat::expect_true(
          SpectralRadius(x = x) > 1
        )
      }
    )
  },
  text = "test-simStateSpace-spectral-radius"
)
