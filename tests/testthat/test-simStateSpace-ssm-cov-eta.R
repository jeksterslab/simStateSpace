## ---- test-simStateSpace-ssm-cov-eta
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
        beta <- matrix(
          data = c(0.7, 0.5, -0.1, 0.0, 0.6, 0.4, 0, 0, 0.5),
          nrow = 3
        )
        psi <- matrix(
          data = c(0.3, 0.1, 0, 0.1, 0.2, 0.05, 0, 0.05, 0.4),
          nrow = 3
        )
        cov_eta <- SSMCovEta(
          beta = beta,
          psi = psi
        )
        testthat::expect_true(
          isSymmetric(cov_eta)
        )
        # cov = beta cov beta' + psi
        testthat::expect_true(
          all(
            abs(cov_eta - (beta %*% cov_eta %*% t(beta) + psi)) <= tol
          )
        )
        testthat::expect_true(
          all(eigen(cov_eta, symmetric = TRUE)$values > 0)
        )
      }
    )
  },
  text = "test-simStateSpace-ssm-cov-eta",
  tol = 1e-10
)
