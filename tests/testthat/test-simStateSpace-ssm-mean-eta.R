## ---- test-simStateSpace-ssm-mean-eta
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
        alpha <- c(0.5, -0.2, 1)
        beta <- matrix(
          data = c(0.7, 0.5, -0.1, 0.0, 0.6, 0.4, 0, 0, 0.5),
          nrow = 3
        )
        mean_eta <- SSMMeanEta(
          beta = beta,
          alpha = alpha
        )
        testthat::expect_true(
          is.null(dim(mean_eta)) && length(mean_eta) == 3
        )
        # mean = alpha + beta mean
        testthat::expect_true(
          all(
            abs(mean_eta - c(alpha + beta %*% mean_eta)) <= tol
          )
        )
        # round trip
        testthat::expect_true(
          all(
            abs(
              SSMInterceptEta(
                beta = beta,
                mean_eta = mean_eta
              ) - alpha
            ) <= tol
          )
        )
      }
    )
  },
  text = "test-simStateSpace-ssm-mean-eta",
  tol = 1e-10
)
