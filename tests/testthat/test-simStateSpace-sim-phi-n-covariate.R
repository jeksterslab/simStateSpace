## ---- test-simStateSpace-sim-phi-n-covariate
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 n) {
    message(text)
    x <- lapply(
      X = seq_len(n),
      FUN = function(i) {
        rnorm(n = 2)
      }
    )
    phi0 <- matrix(
      data = c(
        -0.357, 0.771, -0.450,
        0.0, -0.511, 0.729,
        0, 0, -0.693
      ),
      nrow = 3
    )
    phi1 <- matrix(
      data = 0,
      nrow = 9,
      ncol = 2
    )
    phi1[1, 1] <- 0.01
    phi1[2, 2] <- 0.02
    expected_mean <- colMeans(
      do.call(
        what = "rbind",
        args = lapply(
          X = x,
          FUN = function(x) {
            c(
              c(phi0) + phi1 %*% x
            )
          }
        )
      )
    )
    vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        testthat::expect_true(
          all(
            expected_mean - (1 / n) * Reduce(
              f = `+`,
              x = SimPhiNCovariate(
                n = n,
                phi0 = phi0,
                vcov_phi_vec_l = vcov_phi_vec_l,
                phi1 = phi1,
                x = x
              )
            ) < 0.001
          )
        )
      }
    )
    testthat::test_that(
      paste0(text, "errors"),
      {
        phi_lbound <- matrix(
          data = NA,
          nrow = 2,
          ncol = 2
        )
        diag(phi_lbound) <- -1
        testthat::skip_on_cran()
        testthat::expect_error(
          SimPhiNCovariate(
            n = n,
            phi0 = phi0,
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi1 = phi1,
            x = x,
            phi_lbound = phi_lbound,
            bound = TRUE
          )
        )
        phi_ubound <- matrix(
          data = NA,
          nrow = 2,
          ncol = 2
        )
        diag(phi_ubound) <- 1
        testthat::skip_on_cran()
        testthat::expect_error(
          SimPhiNCovariate(
            n = n,
            phi0 = phi0,
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi1 = phi1,
            x = x,
            phi_ubound = phi_ubound,
            bound = TRUE
          )
        )
        testthat::expect_error(
          SimPhiNCovariate(
            n = n,
            phi0 = diag(3),
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi1 = phi1,
            x = x,
            max_iter = 1
          )
        )
      }
    )
    # coverage
    phi_ubound <- phi_lbound <- matrix(
      data = NA,
      nrow = 3,
      ncol = 3
    )
    diag(phi_lbound) <- -1
    diag(phi_ubound) <- 1
    SimPhiNCovariate(
      n = n,
      phi0 = phi0,
      vcov_phi_vec_l = vcov_phi_vec_l,
      phi1 = phi1,
      x = x,
      phi_lbound = phi_lbound,
      phi_ubound = phi_ubound,
      bound = TRUE
    )
  },
  text = "test-simStateSpace-sim-phi-n-covariate",
  n = 100000
)
