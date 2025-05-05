## ---- test-simStateSpace-ssm-cov
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
        iota <- c(0.317, 0.230)
        phi <- matrix(
          data = c(
            -0.10,
            0.05,
            0.05,
            -0.10
          ),
          nrow = 2
        )
        sigma_l <- t(
          chol(
            matrix(
              data = c(
                2.79,
                0.06,
                0.06,
                3.27
              ),
              nrow = 2
            )
          )
        )
        delta_t <- 1
        ssm <- LinSDE2SSM(
          iota = iota,
          phi = phi,
          sigma_l = sigma_l,
          delta_t = delta_t
        )
        testthat::expect_true(
          all(
            (
              c(
                SSMCovEta(
                  beta = ssm$beta,
                  psi = ssm$psi_l %*% t(ssm$psi_l)
                )
              ) - c(
                19.2,
                10.5,
                10.5,
                21.6
              )
            ) < tol
          )
        )
      }
    )
  },
  text = "test-simStateSpace-ssm-cov",
  tol = 0.000001
)
