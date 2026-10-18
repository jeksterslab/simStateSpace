## ---- test-simStateSpace-lin-sde-cov-y
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        lambda <- matrix(
          data = c(1, 0.8, 0, 0, 0, 0, 1, 0.6, 0, 0, 0, 1),
          nrow = 4
        )
        theta <- diag(c(0.2, 0.3, 0.4, 0.5))
        cov_eta <- matrix(
          data = c(1, 0.3, 0, 0.3, 0.5, 0.1, 0, 0.1, 0.8),
          nrow = 3
        )
        cov_y <- LinSDECovY(
          lambda = lambda,
          theta = theta,
          cov_eta = cov_eta
        )
        testthat::expect_true(
          isSymmetric(cov_y)
        )
        # cov_y = lambda cov_eta lambda' + theta
        testthat::expect_true(
          all(
            abs(cov_y - (lambda %*% cov_eta %*% t(lambda) + theta)) <= tol
          )
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-cov-y",
  tol = 1e-10
)
