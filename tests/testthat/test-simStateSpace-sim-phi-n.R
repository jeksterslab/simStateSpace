## ---- test-simStateSpace-sim-phi-n
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    phi <- matrix(
      data = c(-0.357, 0.771, -0.450, 0.0, -0.511, 0.729, 0, 0, -0.693),
      nrow = 3
    )
    vcov_phi_vec_l <- t(chol(0.01 * diag(9)))
    testthat::test_that(
      paste(text, "hurwitz"),
      {
        testthat::skip_on_cran()
        n <- 20
        for (margin in c(0, -0.3)) {
          set.seed(42)
          draws <- SimPhiN(
            n = n,
            phi = phi,
            vcov_phi_vec_l = vcov_phi_vec_l,
            margin = margin
          )
          set.seed(42)
          z <- .RefStdNormal(n = 50 * n, q = 9)
          ref <- .RefDraw(
            n = n,
            z = z,
            make = function(i, z_i) {
              phi + matrix(data = vcov_phi_vec_l %*% z_i, nrow = 3)
            },
            accept = function(x) {
              .RefTestPhi(x = x, margin = margin)
            }
          )
          if (margin < 0) {
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
        phi_lbound <- matrix(data = NA, nrow = 3, ncol = 3)
        phi_lbound[1, 1] <- -0.357
        phi_ubound <- matrix(data = NA, nrow = 3, ncol = 3)
        phi_ubound[2, 1] <- 0.771
        set.seed(42)
        draws <- SimPhiN(
          n = n,
          phi = phi,
          vcov_phi_vec_l = vcov_phi_vec_l,
          phi_lbound = phi_lbound,
          phi_ubound = phi_ubound,
          bound = TRUE
        )
        set.seed(42)
        z <- .RefStdNormal(n = 50 * n, q = 9)
        ref <- .RefDraw(
          n = n,
          z = z,
          make = function(i, z_i) {
            phi + matrix(data = vcov_phi_vec_l %*% z_i, nrow = 3)
          },
          accept = function(x) {
            .RefInBounds(x = x, lbound = phi_lbound, ubound = phi_ubound) &&
              .RefTestPhi(x = x)
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
          SimPhiN(
            n = 5,
            phi = phi,
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi_lbound = diag(2),
            bound = TRUE
          ),
          regexp = "phi_lbound dims must match phi."
        )
        testthat::expect_error(
          SimPhiN(
            n = 5,
            phi = phi,
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi_ubound = diag(2),
            bound = TRUE
          ),
          regexp = "phi_ubound dims must match phi."
        )
        testthat::expect_error(
          SimPhiN(
            n = 5,
            phi = diag(3),
            vcov_phi_vec_l = vcov_phi_vec_l,
            max_iter = 1
          ),
          regexp = "exceeded max_iter"
        )
      }
    )
  },
  text = "test-simStateSpace-sim-phi-n",
  tol = 1e-12
)
