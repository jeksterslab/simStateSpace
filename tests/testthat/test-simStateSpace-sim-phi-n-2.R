## ---- test-simStateSpace-sim-phi-n-2
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
        phi <- matrix(
          data = c(-0.357, 0.771, -0.450, 0.0, -0.511, 0.729, 0, 0, -0.693),
          nrow = 3
        )
        # large enough for some draws to be projected
        vcov_phi_vec_l <- t(chol(0.04 * diag(9)))
        for (margin in c(1e-3, 0.2)) {
          set.seed(42)
          draws <- SimPhiN2(
            n = n,
            phi = phi,
            vcov_phi_vec_l = vcov_phi_vec_l,
            margin = margin
          )
          set.seed(42)
          z <- .RefStdNormal(n = n, q = 9)
          raw <- .RefDraw(
            n = n,
            z = z,
            make = function(i, z_i) {
              phi + matrix(data = vcov_phi_vec_l %*% z_i, nrow = 3)
            }
          )
          attr(raw, "attempts") <- NULL
          testthat::expect_false(
            all(
              vapply(
                X = raw,
                FUN = .RefHurwitz,
                FUN.VALUE = logical(1),
                margin = margin
              )
            )
          )
          ref <- lapply(
            X = raw,
            FUN = .RefProjectToHurwitz,
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
                FUN = .RefHurwitz,
                FUN.VALUE = logical(1),
                margin = margin
              )
            )
          )
        }
      }
    )
  },
  text = "test-simStateSpace-sim-phi-n-2",
  tol = 1e-10
)
