## ---- test-simStateSpace-log-mat
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      paste(text, "inverse of ExpMat"),
      {
        testthat::skip_on_cran()
        # complex eigenvalues -0.4 +/- 0.6i, inside the principal strip
        phi <- matrix(data = c(-0.4, -0.6, 0.6, -0.4), nrow = 2)
        testthat::expect_true(
          all(abs(LogMat(x = ExpMat(x = phi)) - phi) <= tol)
        )
      }
    )
    testthat::test_that(
      paste(text, "errors"),
      {
        testthat::skip_on_cran()
        testthat::expect_error(
          LogMat(x = matrix(data = 1, nrow = 2, ncol = 3)),
          regexp = "`x` must be a square matrix."
        )
        testthat::expect_error(
          LogMat(x = matrix(data = c(1, 0, 0, Inf), nrow = 2)),
          regexp = "`x` must contain only finite values."
        )
        testthat::expect_error(
          LogMat(x = diag(2), tol = Inf),
          regexp = "`tol` must be finite."
        )
        testthat::expect_error(
          LogMat(x = diag(2), tol = 0),
          regexp = "`tol` must be positive."
        )
        # negative eigenvalue: no real logarithm
        testthat::expect_error(
          LogMat(x = diag(c(-1, 1))),
          regexp = "non-negligible imaginary component"
        )
        testthat::expect_error(
          ExpMat(x = matrix(data = 1, nrow = 2, ncol = 3)),
          regexp = "`x` must be a square matrix."
        )
      }
    )
  },
  text = "test-simStateSpace-log-mat",
  tol = 1e-12
)
