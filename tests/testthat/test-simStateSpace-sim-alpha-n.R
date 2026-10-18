## ---- test-simStateSpace-sim-alpha-n
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
        alpha <- c(1, -0.5, 0.2)
        vcov_alpha_l <- t(
          chol(
            matrix(data = c(1, 0.3, 0, 0.3, 0.5, 0.1, 0, 0.1, 0.8), nrow = 3)
          )
        )
        set.seed(42)
        draws <- SimAlphaN(
          n = n,
          alpha = alpha,
          vcov_alpha_l = vcov_alpha_l
        )
        set.seed(42)
        z <- .RefStdNormal(n = n, q = 3)
        ref <- .RefDraw(
          n = n,
          z = z,
          make = function(i, z_i) {
            c(alpha + vcov_alpha_l %*% z_i)
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
  text = "test-simStateSpace-sim-alpha-n",
  tol = 1e-12
)
