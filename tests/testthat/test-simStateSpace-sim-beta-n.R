## ---- test-simStateSpace-sim-beta-n
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    beta <- matrix(
      data = c(0.7, 0.5, -0.1, 0.0, 0.6, 0.4, 0, 0, 0.5),
      nrow = 3
    )
    vcov_beta_vec_l <- t(chol(0.01 * diag(9)))
    testthat::test_that(
      paste(text, "stationarity"),
      {
        testthat::skip_on_cran()
        n <- 20
        for (margin in c(1, 0.8)) {
          set.seed(42)
          draws <- SimBetaN(
            n = n,
            beta = beta,
            vcov_beta_vec_l = vcov_beta_vec_l,
            margin = margin
          )
          set.seed(42)
          z <- .RefStdNormal(n = 50 * n, q = 9)
          ref <- .RefDraw(
            n = n,
            z = z,
            make = function(i, z_i) {
              beta + matrix(data = vcov_beta_vec_l %*% z_i, nrow = 3)
            },
            accept = function(x) {
              .RefStationary(x = x, margin = margin)
            }
          )
          if (margin < 1) {
            # the rejection step was used
            testthat::expect_true(
              attr(ref, "attempts") > n
            )
          }
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
      paste(text, "bounds"),
      {
        testthat::skip_on_cran()
        n <- 20
        beta_lbound <- matrix(data = NA, nrow = 3, ncol = 3)
        beta_lbound[1, 1] <- 0.7
        beta_ubound <- matrix(data = NA, nrow = 3, ncol = 3)
        beta_ubound[2, 1] <- 0.5
        set.seed(42)
        draws <- SimBetaN(
          n = n,
          beta = beta,
          vcov_beta_vec_l = vcov_beta_vec_l,
          beta_lbound = beta_lbound,
          beta_ubound = beta_ubound,
          bound = TRUE
        )
        set.seed(42)
        z <- .RefStdNormal(n = 50 * n, q = 9)
        ref <- .RefDraw(
          n = n,
          z = z,
          make = function(i, z_i) {
            beta + matrix(data = vcov_beta_vec_l %*% z_i, nrow = 3)
          },
          accept = function(x) {
            .RefInBounds(x = x, lbound = beta_lbound, ubound = beta_ubound) &&
              .RefStationary(x = x)
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
    )
    testthat::test_that(
      paste(text, "errors"),
      {
        testthat::skip_on_cran()
        testthat::expect_error(
          SimBetaN(
            n = 5,
            beta = beta,
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta_lbound = diag(2),
            bound = TRUE
          ),
          regexp = "beta_lbound dims must match beta."
        )
        testthat::expect_error(
          SimBetaN(
            n = 5,
            beta = beta,
            vcov_beta_vec_l = vcov_beta_vec_l,
            beta_ubound = diag(2),
            bound = TRUE
          ),
          regexp = "beta_ubound dims must match beta."
        )
        testthat::expect_error(
          SimBetaN(
            n = 5,
            beta = 5 * diag(3),
            vcov_beta_vec_l = vcov_beta_vec_l,
            max_iter = 1
          ),
          regexp = "exceeded max_iter"
        )
      }
    )
  },
  text = "test-simStateSpace-sim-beta-n",
  tol = 1e-12
)
