## ---- test-simStateSpace-sim-cov-n
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
        sigma <- matrix(
          data = c(1.0, 0.5, 0.3, 0.5, 1.0, 0.4, 0.3, 0.4, 1.0),
          nrow = 3
        )
        # large enough for some draws to be clipped to positive definite
        vcov_sigma_vech_l <- t(chol(0.25 * diag(6)))
        set.seed(42)
        draws <- SimCovN(
          n = n,
          sigma = sigma,
          vcov_sigma_vech_l = vcov_sigma_vech_l
        )
        set.seed(42)
        z <- .RefStdNormal(n = n, q = 6)
        clipped <- 0
        ref <- .RefDraw(
          n = n,
          z = z,
          make = function(i, z_i) {
            vech_i <- sigma[lower.tri(sigma, diag = TRUE)] +
              c(vcov_sigma_vech_l %*% z_i)
            sigma_i <- matrix(data = 0, nrow = 3, ncol = 3)
            sigma_i[lower.tri(sigma_i, diag = TRUE)] <- vech_i
            sigma_i[upper.tri(sigma_i)] <- t(sigma_i)[upper.tri(sigma_i)]
            eig <- eigen(sigma_i, symmetric = TRUE)
            if (any(eig$values < 1e-8)) {
              clipped <<- clipped + 1
            }
            eig$vectors %*% diag(pmax(eig$values, 1e-8)) %*% t(eig$vectors)
          }
        )
        attr(ref, "attempts") <- NULL
        testthat::expect_true(
          clipped > 0
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
              FUN = function(x) {
                isSymmetric(x) &&
                  min(eigen(x, symmetric = TRUE)$values) > 0
              },
              FUN.VALUE = logical(1)
            )
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-cov-n",
  tol = 1e-10
)
