## ---- test-simStateSpace-spectral-abscissa
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      paste0(text, "x < 0"),
      {
        testthat::skip_on_cran()
        x <- matrix(
          data = c(
            -0.5, -0.2,
            1.0, -0.3
          ),
          nrow = 2
        )
        testthat::expect_true(
          SpectralAbscissa(x = x) < 0
        )
      }
    )
    testthat::test_that(
      paste0(text, "x > 0"),
      {
        testthat::skip_on_cran()
        x <- matrix(
          data = c(
            0.10, 0.50,
            -0.40, 0.20
          ),
          nrow = 2
        )
        testthat::expect_true(
          SpectralAbscissa(x = x) > 0
        )
      }
    )
  },
  text = "test-simStateSpace-spectral-abscissa"
)
