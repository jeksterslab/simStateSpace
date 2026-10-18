## ---- test-simStateSpace-sim-phi-n-covariate
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    n <- 20
    phi0 <- matrix(
      data = c(-0.357, 0.771, -0.450, 0.0, -0.511, 0.729, 0, 0, -0.693),
      nrow = 3
    )
    phi1 <- matrix(data = 0, nrow = 9, ncol = 2)
    phi1[1, 1] <- 0.1
    phi1[2, 2] <- -0.2
    phi1[5, 1] <- 0.05
    vcov_phi_vec_l <- t(chol(0.01 * diag(9)))
    testthat::test_that(
      paste(text, "hurwitz and bounds"),
      {
        testthat::skip_on_cran()
        set.seed(1)
        x <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            stats::rnorm(n = 2)
          }
        )
        phi_lbound <- matrix(data = NA, nrow = 3, ncol = 3)
        phi_lbound[1, 1] <- -0.4
        phi_ubound <- matrix(data = NA, nrow = 3, ncol = 3)
        phi_ubound[2, 1] <- 1.2
        for (bound in c(FALSE, TRUE)) {
          set.seed(42)
          draws <- SimPhiNCovariate(
            n = n,
            phi0 = phi0,
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi1 = phi1,
            x = x,
            margin = -0.3,
            phi_lbound = phi_lbound,
            phi_ubound = phi_ubound,
            bound = bound
          )
          set.seed(42)
          z <- .RefStdNormal(n = 100 * n, q = 9)
          ref <- .RefDraw(
            n = n,
            z = z,
            make = function(i, z_i) {
              matrix(
                data = c(phi0) + phi1 %*% x[[i]] + vcov_phi_vec_l %*% z_i,
                nrow = 3
              )
            },
            accept = function(x) {
              (!bound || .RefInBounds(
                x = x,
                lbound = phi_lbound,
                ubound = phi_ubound
              )) &&
                .RefTestPhi(x = x, margin = -0.3)
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
          SimPhiNCovariate(
            n = 5,
            phi0 = phi0,
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi1 = phi1,
            x = x,
            phi_lbound = diag(2),
            bound = TRUE
          ),
          regexp = "phi_lbound dims must match phi0."
        )
        testthat::expect_error(
          SimPhiNCovariate(
            n = 5,
            phi0 = phi0,
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi1 = phi1,
            x = x,
            phi_ubound = diag(2),
            bound = TRUE
          ),
          regexp = "phi_ubound dims must match phi0."
        )
        testthat::expect_error(
          SimPhiNCovariate(
            n = 5,
            phi0 = diag(3),
            vcov_phi_vec_l = vcov_phi_vec_l,
            phi1 = phi1,
            x = x,
            max_iter = 1
          ),
          regexp = "exceeded max_iter"
        )
      }
    )
  },
  text = "test-simStateSpace-sim-phi-n-covariate",
  tol = 1e-12
)
