## ---- test-simStateSpace-ssm-2-lin-sde-edge
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      paste(text, "rank-deficient covariance"),
      {
        testthat::skip_on_cran()
        # rank-2 diffusion in three dimensions: chol() fails on psi, and the
        # unpivoted fallback in 1.2.16 rejected it as not semidefinite
        sigma_l <- matrix(
          data = c(
            0.06, -2.42, -0.34,
            0, 0.03, -0.39,
            0, 0, 0
          ),
          nrow = 3
        )
        sigma <- tcrossprod(sigma_l)
        phi <- diag(x = -0.5, nrow = 3)
        ssm <- simStateSpace::LinSDE2SSM(
          iota = rep(x = 0, times = 3),
          phi = phi,
          sigma_l = sigma_l,
          delta_t = 1
        )
        testthat::expect_true(
          all(ssm$psi_l[upper.tri(ssm$psi_l)] == 0)
        )
        testthat::expect_true(
          all(diag(ssm$psi_l) >= 0)
        )
        testthat::expect_true(
          max(
            abs(
              tcrossprod(ssm$psi_l) - (1 - exp(-1)) * sigma
            )
          ) <= tol
        )
        linsde <- simStateSpace::SSM2LinSDE(
          alpha = ssm$alpha,
          beta = ssm$beta,
          psi_l = ssm$psi_l,
          delta_t = 1
        )
        testthat::expect_true(
          max(
            abs(
              tcrossprod(linsde$sigma_l) - sigma
            )
          ) <= tol
        )
        testthat::expect_true(
          max(
            abs(
              linsde$phi - phi
            )
          ) <= tol
        )
      }
    )
    testthat::test_that(
      paste(text, "singular beta"),
      {
        testthat::skip_on_cran()
        # diagonalizable: the logarithm has -Inf entries
        testthat::expect_error(
          object = simStateSpace::SSM2LinSDE(
            alpha = c(0, 0),
            beta = diag(x = c(0.5, 0)),
            psi_l = diag(x = 2),
            delta_t = 1
          ),
          regexp = "not finite"
        )
        # not diagonalizable: no logarithm
        testthat::expect_error(
          object = simStateSpace::SSM2LinSDE(
            alpha = c(0, 0),
            beta = matrix(data = c(0.5, 0, 1, 0), nrow = 2),
            psi_l = diag(x = 2),
            delta_t = 1
          ),
          regexp = "logarithm of beta failed"
        )
      }
    )
    testthat::test_that(
      paste(text, "tol_psi_l"),
      {
        testthat::skip_on_cran()
        psi_l <- diag(x = 2)
        psi_l[1, 2] <- 1e-06
        testthat::expect_error(
          object = simStateSpace::SSM2LinSDE(
            alpha = c(0, 0),
            beta = diag(x = 0.5, nrow = 2),
            psi_l = psi_l,
            delta_t = 1
          ),
          regexp = "lower triangular"
        )
        # the triangularity tolerance does not loosen the imaginary-part check
        linsde <- simStateSpace::SSM2LinSDE(
          alpha = c(0, 0),
          beta = diag(x = 0.5, nrow = 2),
          psi_l = psi_l,
          delta_t = 1,
          tol_psi_l = 1e-04
        )
        testthat::expect_true(
          max(
            abs(
              linsde$phi - diag(x = log(0.5), nrow = 2)
            )
          ) <= tol
        )
        testthat::expect_error(
          object = simStateSpace::SSM2LinSDE(
            alpha = 0,
            beta = matrix(data = -1),
            psi_l = matrix(data = 1),
            delta_t = 1,
            tol_psi_l = 1
          ),
          regexp = "imaginary"
        )
        testthat::expect_error(
          object = simStateSpace::SSM2LinSDE(
            alpha = 0,
            beta = matrix(data = 0.5),
            psi_l = matrix(data = 1),
            delta_t = 1,
            tol_psi_l = -1
          ),
          regexp = "tol_psi_l"
        )
      }
    )
  },
  text = "test-simStateSpace-ssm-2-lin-sde-edge",
  tol = 1e-12
)
