## ---- test-simStateSpace-lin-sde-mean-y
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
        nu <- c(0.5, 1, -0.5, 0.2)
        lambda <- matrix(
          data = c(1, 0.8, 0, 0, 0, 0, 1, 0.6, 0, 0, 0, 1),
          nrow = 4
        )
        mean_eta <- c(1, -2, 0.5)
        mean_y <- LinSDEMeanY(
          nu = nu,
          lambda = lambda,
          mean_eta = mean_eta
        )
        testthat::expect_true(
          is.null(dim(mean_y)) && length(mean_y) == 4
        )
        # mean_y = nu + lambda mean_eta
        testthat::expect_true(
          all(
            abs(mean_y - c(nu + lambda %*% mean_eta)) <= tol
          )
        )
        # round trip
        testthat::expect_true(
          all(
            abs(
              LinSDEInterceptY(
                mean_y = mean_y,
                mean_eta = mean_eta,
                lambda = lambda
              ) - nu
            ) <= tol
          )
        )
        # column matrices are accepted as vectors
        testthat::expect_identical(
          object = LinSDEMeanY(
            nu = matrix(data = nu),
            lambda = lambda,
            mean_eta = matrix(data = mean_eta)
          ),
          expected = mean_y
        )
        testthat::expect_identical(
          object = LinSDEInterceptY(
            mean_y = matrix(data = mean_y),
            mean_eta = matrix(data = mean_eta),
            lambda = lambda
          ),
          expected = LinSDEInterceptY(
            mean_y = mean_y,
            mean_eta = mean_eta,
            lambda = lambda
          )
        )
      }
    )
  },
  text = "test-simStateSpace-lin-sde-mean-y",
  tol = 1e-10
)
