## ---- test-simStateSpace-sim-cov-diag-n
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
        sigma_diag <- c(1, 0.5, 0.2)
        # large enough for some draws to be clipped at 1e-8
        vcov_sigma_diag_l <- t(chol(0.1 * diag(3)))
        set.seed(42)
        draws <- SimCovDiagN(
          n = n,
          sigma_diag = sigma_diag,
          vcov_sigma_diag_l = vcov_sigma_diag_l
        )
        set.seed(42)
        z <- .RefStdNormal(n = n, q = 3)
        raw <- .RefDraw(
          n = n,
          z = z,
          make = function(i, z_i) {
            sigma_diag + c(vcov_sigma_diag_l %*% z_i)
          }
        )
        testthat::expect_true(
          any(unlist(raw) < 1e-8)
        )
        ref <- lapply(
          X = raw,
          FUN = function(x) {
            diag(pmax(x, 1e-8))
          }
        )
        testthat::expect_equal(
          object = draws,
          expected = ref,
          tolerance = tol
        )
      }
    )
  },
  text = "test-simStateSpace-sim-cov-diag-n",
  tol = 1e-12
)
