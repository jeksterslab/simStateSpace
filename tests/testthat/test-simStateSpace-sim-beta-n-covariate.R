## ---- test-simStateSpace-sim-beta-n-covariate
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
    beta0 <- matrix(
      data = c(
        0.7, 0.5, -0.1,
        0.0, 0.6, 0.4,
        0, 0, 0.5
      ),
      nrow = 3
    )
    beta1 <- matrix(
      data = 0,
      nrow = 9,
      ncol = 2
    )
    beta1[1, 1] <- 0.01
    beta1[2, 2] <- 0.02
    expected_mean <- colMeans(
      do.call(
        what = "rbind",
        args = lapply(
          X = x,
          FUN = function(x) {
            c(
              c(beta0) + beta1 %*% x
            )
          }
        )
      )
    )
    vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        testthat::expect_true(
          all(
            expected_mean - (1 / n) * Reduce(
              f = `+`,
              x = SimBetaNCovariate(
                n = n,
                beta0 = beta0,
                vcov_beta_vec_l = vcov_beta_vec_l,
                beta1 = beta1,
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
        beta_lbound <- matrix(
          data = NA,
          nrow = 2,
          ncol = 2
        )
        diag(beta_lbound) <- -1
        testthat::skip_on_cran()
        testthat::expect_error(
          SimBetaNCovariate(
            n = n,
            beta0 = beta0,
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta1 = beta1,
            x = x,
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
          SimBetaNCovariate(
            n = n,
            beta0 = beta0,
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta1 = beta1,
            x = x,
            beta_ubound = beta_ubound,
            bound = TRUE
          )
        )
        testthat::expect_error(
          SimBetaNCovariate(
            n = n,
            beta0 = 5 * diag(3),
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta1 = beta1,
            x = x,
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
    SimBetaNCovariate(
      n = n,
      beta0 = beta0,
      vcov_beta_vec_l = vcov_beta_vec_l,
      beta1 = beta1,
      x = x,
      beta_lbound = beta_lbound,
      beta_ubound = beta_ubound,
      bound = TRUE
    )
  },
  text = "test-simStateSpace-sim-beta-n-covariate",
  n = 100000
)
