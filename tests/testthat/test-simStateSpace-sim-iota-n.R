## ---- test-simStateSpace-sim-iota-n
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
        iota <- c(1, -0.5, 0.2)
        vcov_iota_l <- t(
          chol(
            matrix(data = c(1, 0.3, 0, 0.3, 0.5, 0.1, 0, 0.1, 0.8), nrow = 3)
          )
        )
        set.seed(42)
        draws <- SimIotaN(
          n = n,
          iota = iota,
          vcov_iota_l = vcov_iota_l
        )
        set.seed(42)
        z <- .RefStdNormal(n = n, q = 3)
        ref <- .RefDraw(
          n = n,
          z = z,
          make = function(i, z_i) {
            c(iota + vcov_iota_l %*% z_i)
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
  text = "test-simStateSpace-sim-iota-n",
  tol = 1e-12
)
