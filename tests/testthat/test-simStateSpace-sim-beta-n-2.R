## ---- test-simStateSpace-sim-beta-n-2
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        n <- 50
        beta <- matrix(
          data = c(0.7, 0.5, -0.1, 0.0, 0.6, 0.4, 0, 0, 0.5),
          nrow = 3
        )
        # large enough for some draws to be projected
        vcov_beta_vec_l <- t(chol(0.04 * diag(9)))
        for (margin in c(0.98, 0.9)) {
          set.seed(42)
          draws <- SimBetaN2(
            n = n,
            beta = beta,
            vcov_beta_vec_l = vcov_beta_vec_l,
            margin = margin
          )
          set.seed(42)
          z <- .RefStdNormal(n = n, q = 9)
          raw <- .RefDraw(
            n = n,
            z = z,
            make = function(i, z_i) {
              beta + matrix(data = vcov_beta_vec_l %*% z_i, nrow = 3)
            }
          )
          attr(raw, "attempts") <- NULL
          testthat::expect_false(
            all(
              vapply(
                X = raw,
                FUN = .RefStationary,
                FUN.VALUE = logical(1),
                margin = margin
              )
            )
          )
          ref <- lapply(
            X = raw,
            FUN = .RefProjectToStability,
            margin = margin
          )
          testthat::expect_equal(
            object = draws,
            expected = ref,
            tolerance = tol
          )
          testthat::expect_true(
            all(
              vapply(
                X = draws,
                FUN = .RefStationary,
                FUN.VALUE = logical(1),
                margin = margin
              )
            )
          )
        }
      }
    )
  },
  text = "test-simStateSpace-sim-beta-n-2",
  tol = 1e-10
)
