## ---- test-simStateSpace-sim-nu-n
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
        n <- 20
        nu <- c(1, -0.5, 0.2)
        vcov_nu_l <- t(
          chol(
            matrix(data = c(1, 0.3, 0, 0.3, 0.5, 0.1, 0, 0.1, 0.8), nrow = 3)
          )
        )
        set.seed(42)
        draws <- SimNuN(
          n = n,
          nu = nu,
          vcov_nu_l = vcov_nu_l
        )
        set.seed(42)
        z <- .RefStdNormal(n = n, q = 3)
        ref <- .RefDraw(
          n = n,
          z = z,
          make = function(i, z_i) {
            c(nu + vcov_nu_l %*% z_i)
          }
        )
        attr(ref, "attempts") <- NULL
        testthat::expect_equal(
          object = draws,
          expected = ref,
          tolerance = tol
        )
      }
    )
  },
  text = "test-simStateSpace-sim-nu-n",
  tol = 1e-12
)
