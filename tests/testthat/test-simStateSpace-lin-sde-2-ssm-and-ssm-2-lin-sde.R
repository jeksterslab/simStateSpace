## ---- test-simStateSpace-lin-sde-2-ssm-and-ssm-2-lin-sde
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol,
                 tol_chow) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        delta_t <- 0.2
        k <- p <- 2
        iota <- c(
          0.32,
          0.23
        )
        phi <- matrix(
          data = c(
            -0.10,
            0.05,
            0.05,
            -0.10
          ),
          nrow = k
        )
        sigma <- matrix(
          data = c(
            2.79,
            0.06,
            0.06,
            3.27
          ),
          nrow = p
        )
        sigma_l <- t(chol(sigma))
        I <- diag(p)
        beta <- simStateSpace::ExpMat(
          x = phi * delta_t
        )
        phi_log <- simStateSpace::LogMat(
          x = beta
        ) / delta_t
        alpha <- solve(
          a = phi,
          b = (beta - I) %*% iota
        )
        phi_hashtag <- kronecker(phi, I) + kronecker(I, phi)
        J <- diag(p * p)
        psi_vec <- solve(
          a = phi_hashtag,
          b = (
            simStateSpace::ExpMat(
              x = phi_hashtag * delta_t
            ) - J
          ) %*% c(sigma)
        )
        psi <- matrix(
          data = psi_vec,
          nrow = p
        )
        psi <- (psi + t(psi)) / 2
        psi_l <- t(chol(psi))
        ssm <- simStateSpace::LinSDE2SSM(
          iota = iota,
          phi = phi,
          sigma_l = sigma_l,
          delta_t = delta_t
        )
        testthat::expect_true(
          all(
            abs(
              c(phi_log) - c(phi)
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              c(ssm$alpha) - c(alpha)
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              c(ssm$beta) - c(beta)
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              c(ssm$psi_l) - c(psi_l)
            ) <= tol
          )
        )
        linsde <- simStateSpace::SSM2LinSDE(
          alpha = alpha,
          beta = beta,
          psi_l = psi_l,
          delta_t = delta_t
        )
        testthat::expect_true(
          all(
            abs(
              c(linsde$iota) - c(iota)
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              c(linsde$phi) - c(phi)
            ) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              c(linsde$sigma_l) - c(sigma_l)
            ) <= tol
          )
        )
        # Rounded values for delta_t = 0.2.
        # Chow et al. Table 2 reports beta and psi.
        alpha_chow <- c(
          0.064,
          0.046
        )
        beta_chow <- matrix(
          data = c(
            0.98,
            0.01,
            0.01,
            0.98
          ),
          nrow = p
        )
        psi_chow <- matrix(
          data = c(
            0.55,
            0.02,
            0.02,
            0.64
          ),
          nrow = p
        )
        testthat::expect_true(
          all(
            abs(
              c(ssm$alpha) - c(alpha_chow)
            ) <= tol_chow
          )
        )
        testthat::expect_true(
          all(
            abs(
              c(ssm$beta) - c(beta_chow)
            ) <= tol_chow
          )
        )
        testthat::expect_true(
          all(
            abs(
              c(psi) - c(psi_chow)
            ) <= tol_chow
          )
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-2-ssm-and-ssm-2-lin-sde",
  tol = 1.0e-08,
  tol_chow = 1.0e-02
)
