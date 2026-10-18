## ---- test-simStateSpace-lin-sde-mean-eta
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
        iota <- c(1, -0.5, 0.3)
        phi <- matrix(
          data = c(-0.357, 0.771, -0.450, 0.0, -0.511, 0.729, 0, 0, -0.693),
          nrow = 3
        )
        mean_eta <- LinSDEMeanEta(
          phi = phi,
          iota = iota
        )
        testthat::expect_true(
          is.null(dim(mean_eta)) && length(mean_eta) == 3
        )
        # phi mean + iota = 0
        testthat::expect_true(
          all(
            abs(phi %*% mean_eta + iota) <= tol
          )
        )
        # round trip
        testthat::expect_true(
          all(
            abs(
              LinSDEInterceptEta(
                phi = phi,
                mean_eta = mean_eta
              ) - iota
            ) <= tol
          )
        )
        # column matrices are accepted as vectors
        testthat::expect_identical(
          object = LinSDEMeanEta(
            phi = phi,
            iota = matrix(data = iota)
          ),
          expected = mean_eta
        )
        testthat::expect_identical(
          object = LinSDEInterceptEta(
            phi = phi,
            mean_eta = matrix(data = mean_eta)
          ),
          expected = LinSDEInterceptEta(
            phi = phi,
            mean_eta = mean_eta
          )
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-mean-eta",
  tol = 1e-10
)
