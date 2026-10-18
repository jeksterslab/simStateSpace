## ---- test-simStateSpace-lin-sde-cov-delta-t
lapply(
  X = c(
    0.2,
    1.0,
    5.0
  ),
  FUN = function(delta_t,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      paste(text, delta_t),
      {
        testthat::skip_on_cran()
        # scalar OU: closed form
        testthat::expect_equal(
          object = c(
            simStateSpace::LinSDECovDeltaT(
              phi = matrix(data = -0.4),
              sigma = matrix(data = 1.7),
              delta_t = delta_t
            )
          ),
          expected = 1.7 * expm1(x = -0.8 * delta_t) / -0.8,
          tolerance = tol
        )
        # zero drift: delta_t * sigma
        sigma <- matrix(
          data = c(
            1,
            0.2,
            0.2,
            0.6
          ),
          nrow = 2
        )
        testthat::expect_equal(
          object = simStateSpace::LinSDECovDeltaT(
            phi = matrix(
              data = 0,
              nrow = 2,
              ncol = 2
            ),
            sigma = sigma,
            delta_t = delta_t
          ),
          expected = delta_t * sigma,
          tolerance = tol
        )
        # agrees with the process noise of LinSDE2SSM
        phi <- matrix(
          data = c(
            -0.3,
            0,
            0.7,
            -0.8
          ),
          nrow = 2
        )
        ssm <- simStateSpace::LinSDE2SSM(
          iota = c(
            0,
            0
          ),
          phi = phi,
          sigma_l = t(x = chol(x = sigma)),
          delta_t = delta_t
        )
        psi <- simStateSpace::LinSDECovDeltaT(
          phi = phi,
          sigma = sigma,
          delta_t = delta_t
        )
        testthat::expect_equal(
          object = psi,
          expected = tcrossprod(x = ssm$psi_l),
          tolerance = tol
        )
        # stationary covariance: Cov(eta) = beta Cov(eta) beta' + psi
        cov_eta <- simStateSpace::LinSDECovEta(
          phi = phi,
          sigma = sigma
        )
        testthat::expect_equal(
          object = psi,
          expected = cov_eta - ssm$beta %*% cov_eta %*% t(x = ssm$beta),
          tolerance = tol
        )
        # linear in sigma, so an indefinite sigma is allowed
        basis <- matrix(
          data = c(
            0,
            1,
            1,
            0
          ),
          nrow = 2
        )
        testthat::expect_equal(
          object = simStateSpace::LinSDECovDeltaT(
            phi = phi,
            sigma = sigma - basis,
            delta_t = delta_t
          ),
          expected = psi - simStateSpace::LinSDECovDeltaT(
            phi = phi,
            sigma = basis,
            delta_t = delta_t
          ),
          tolerance = tol
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-cov-delta-t",
  tol = 1.0e-10
)

lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        testthat::expect_error(
          object = simStateSpace::LinSDECovDeltaT(
            phi = matrix(data = -0.4),
            sigma = matrix(
              data = 1,
              nrow = 2,
              ncol = 2
            )
          ),
          regexp = "same dimensions"
        )
        testthat::expect_error(
          object = simStateSpace::LinSDECovDeltaT(
            phi = diag(x = -0.4, nrow = 2),
            sigma = matrix(
              data = c(
                1,
                0.5,
                0,
                1
              ),
              nrow = 2
            )
          ),
          regexp = "symmetric"
        )
        testthat::expect_error(
          object = simStateSpace::LinSDECovDeltaT(
            phi = matrix(data = -0.4),
            sigma = matrix(data = 1),
            delta_t = 0
          ),
          regexp = "positive"
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-cov-delta-t-errors"
)
