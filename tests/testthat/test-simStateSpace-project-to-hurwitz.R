## ---- test-simStateSpace-project-to-hurwitz
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      paste0(text, "< 0"),
      {
        testthat::skip_on_cran()
        x <- matrix(
          data = c(
            -0.50, -0.20,
            1.00, -0.30
          ),
          nrow = 2
        )
        testthat::expect_true(
          SpectralAbscissa(x = x) < 0
        )
        testthat::expect_true(
          identical(ProjectToHurwitz(x = x), x)
        )
      }
    )
    testthat::test_that(
      paste0(text, ">= 0"),
      {
        testthat::skip_on_cran()
        x <- matrix(
          data = c(
            0.10, -0.40,
            0.50, 0.20
          ),
          nrow = 2
        )
        testthat::expect_true(
          SpectralAbscissa(x = x) >= 0
        )
        testthat::expect_true(
          SpectralAbscissa(x = ProjectToHurwitz(x = x)) <= -1e-3
        )
      }
    )
  },
  text = "test-simStateSpace-project-to-hurwitz"
)
