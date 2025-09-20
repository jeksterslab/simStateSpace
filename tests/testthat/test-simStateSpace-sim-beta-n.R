## ---- test-simStateSpace-sim-beta-n
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 n) {
    message(text)
    beta <- matrix(
      data = c(
        0.7, 0.5, -0.1,
        0.0, 0.6, 0.4,
        0, 0, 0.5
      ),
      nrow = 3
    )
    vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        testthat::expect_true(
          all(
            beta - (1 / n) * Reduce(
              f = `+`,
              x = SimBetaN(
                n = n,
                beta = beta,
                vcov_beta_vec_l = vcov_beta_vec_l
              )
            ) < 0.001
          )
        )
      }
    )
    testthat::test_that(
      paste0(text, "errors"),
      {
        beta_lbound <- matrix(
          data = NA,
          nrow = 2,
          ncol = 2
        )
        diag(beta_lbound) <- -1
        testthat::skip_on_cran()
        testthat::expect_error(
          SimBetaN(
            n = n,
            beta = beta,
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta_lbound = beta_lbound,
            bound = TRUE
          )
        )
        beta_ubound <- matrix(
          data = NA,
          nrow = 2,
          ncol = 2
        )
        diag(beta_ubound) <- 1
        testthat::skip_on_cran()
        testthat::expect_error(
          SimBetaN(
            n = n,
            beta = beta,
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta_ubound = beta_ubound,
            bound = TRUE
          )
        )
        testthat::expect_error(
          SimBetaN(
            n = n,
            beta = 5 * diag(3),
            vcov_beta_vec_l = vcov_beta_vec_l,
            max_iter = 1
          )
        )
      }
    )
    # coverage
    beta_ubound <- beta_lbound <- matrix(
      data = NA,
      nrow = 3,
      ncol = 3
    )
    diag(beta_lbound) <- -1
    diag(beta_ubound) <- 1
    SimBetaN(
      n = n,
      beta = beta,
      vcov_beta_vec_l = vcov_beta_vec_l,
      beta_lbound = beta_lbound,
      beta_ubound = beta_ubound,
      bound = TRUE
    )
  },
  text = "test-simStateSpace-sim-beta-n",
  n = 100000
)
