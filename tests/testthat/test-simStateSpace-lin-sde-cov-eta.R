## ---- test-simStateSpace-lin-sde-cov-eta
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
        phi <- matrix(
          data = c(-0.357, 0.771, -0.450, 0.0, -0.511, 0.729, 0, 0, -0.693),
          nrow = 3
        )
        sigma <- matrix(
          data = c(
            0.24455556, 0.02201587, -0.05004762,
            0.02201587, 0.07067800, 0.01539456,
            -0.05004762, 0.01539456, 0.07553061
          ),
          nrow = 3
        )
        cov_eta <- LinSDECovEta(
          phi = phi,
          sigma = sigma
        )
        testthat::expect_true(
          isSymmetric(cov_eta)
        )
        # phi cov + cov phi' + sigma = 0
        testthat::expect_true(
          all(
            abs(phi %*% cov_eta + cov_eta %*% t(phi) + sigma) <= tol
          )
        )
        # internal solvers: A X + X B + C = 0
        syl <- .SolveSyl(
          A = phi,
          B = t(phi),
          C = sigma
        )
        testthat::expect_true(
          all(
            abs(phi %*% syl + syl %*% t(phi) + sigma) <= tol
          )
        )
        testthat::expect_true(
          all(
            abs(
              .SolveLya(
                A = phi,
                Q = sigma
              ) - cov_eta
            ) <= tol
          )
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-cov-eta",
  tol = 1e-10
)
