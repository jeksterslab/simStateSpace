## ---- test-simStateSpace-sim-phi-n
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 n) {
    message(text)
    phi <- matrix(
      data = c(
        -0.357, 0.771, -0.450,
        0.0, -0.511, 0.729,
        0, 0, -0.693
      ),
      nrow = 3
    )
    vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        testthat::expect_true(
          all(
            phi - (1 / n) * Reduce(
              f = `+`,
              x = SimPhiN(
                n = n,
                phi = phi,
                vcov_phi_vec_l = vcov_phi_vec_l
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
          SimPhiN(
            n = n,
            phi = phi,
            vcov_phi_vec_l = vcov_phi_vec_l,
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
          SimPhiN(
            n = n,
            phi = phi,
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi_ubound = phi_ubound,
            bound = TRUE
          )
        )
        testthat::expect_error(
          SimPhiN(
            n = n,
            phi = diag(3),
            vcov_phi_vec_l = vcov_phi_vec_l,
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
    SimPhiN(
      n = n,
      phi = phi,
      vcov_phi_vec_l = vcov_phi_vec_l,
      phi_lbound = phi_lbound,
      phi_ubound = phi_ubound,
      bound = TRUE
    )
  },
  text = "test-simStateSpace-sim-phi-n",
  n = 100000
)
