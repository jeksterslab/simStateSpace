## ---- test-simStateSpace-sim-beta-n-covariate
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    n <- 20
    beta0 <- matrix(
      data = c(0.7, 0.5, -0.1, 0.0, 0.6, 0.4, 0, 0, 0.5),
      nrow = 3
    )
    beta1 <- matrix(data = 0, nrow = 9, ncol = 2)
    beta1[1, 1] <- 0.1
    beta1[2, 2] <- -0.2
    beta1[5, 1] <- 0.05
    vcov_beta_vec_l <- t(chol(0.01 * diag(9)))
    testthat::test_that(
      paste(text, "stationarity and bounds"),
      {
        testthat::skip_on_cran()
        set.seed(1)
        x <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            stats::rnorm(n = 2)
          }
        )
        beta_lbound <- matrix(data = NA, nrow = 3, ncol = 3)
        beta_lbound[1, 1] <- 0.6
        beta_ubound <- matrix(data = NA, nrow = 3, ncol = 3)
        beta_ubound[2, 1] <- 1
        for (bound in c(FALSE, TRUE)) {
          set.seed(42)
          draws <- SimBetaNCovariate(
            n = n,
            beta0 = beta0,
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta1 = beta1,
            x = x,
            margin = 0.9,
            beta_lbound = beta_lbound,
            beta_ubound = beta_ubound,
            bound = bound
          )
          set.seed(42)
          z <- .RefStdNormal(n = 50 * n, q = 9)
          ref <- .RefDraw(
            n = n,
            z = z,
            make = function(i, z_i) {
              matrix(
                data = c(beta0) + beta1 %*% x[[i]] + vcov_beta_vec_l %*% z_i,
                nrow = 3
              )
            },
            accept = function(x) {
              (!bound || .RefInBounds(
                x = x,
                lbound = beta_lbound,
                ubound = beta_ubound
              )) &&
                .RefStationary(x = x, margin = 0.9)
            }
          )
          testthat::expect_true(
            attr(ref, "attempts") > n
          )
          attr(ref, "attempts") <- NULL
          testthat::expect_equal(
            object = draws,
            expected = ref,
            tolerance = tol
          )
        }
      }
    )
    testthat::test_that(
      paste(text, "errors"),
      {
        testthat::skip_on_cran()
        x <- rep(x = list(c(0, 0)), times = 5)
        testthat::expect_error(
          SimBetaNCovariate(
            n = 5,
            beta0 = beta0,
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta1 = beta1,
            x = x,
            beta_lbound = diag(2),
            bound = TRUE
          ),
          regexp = "beta_lbound dims must match beta0."
        )
        testthat::expect_error(
          SimBetaNCovariate(
            n = 5,
            beta0 = beta0,
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta1 = beta1,
            x = x,
            beta_ubound = diag(2),
            bound = TRUE
          ),
          regexp = "beta_ubound dims must match beta0."
        )
        testthat::expect_error(
          SimBetaNCovariate(
            n = 5,
            beta0 = 5 * diag(3),
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta1 = beta1,
            x = x,
            max_iter = 1
          ),
          regexp = "exceeded max_iter"
        )
      }
    )
  },
  text = "test-simStateSpace-sim-beta-n-covariate",
  tol = 1e-12
)
