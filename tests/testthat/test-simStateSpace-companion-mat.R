## ---- test-simStateSpace-companion-mat
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      "CompanionMat() returns the input matrix for order 1.",
      {
        testthat::skip_on_cran()
        beta <- list(
          matrix(
            data = c(
              0.40, 0.10,
              0.20, 0.30
            ),
            nrow = 2,
            ncol = 2
          )
        )

        out <- CompanionMat(x = beta)

        testthat::expect_equal(
          object = out,
          expected = beta[[1]]
        )
      }
    )

    testthat::test_that(
      "CompanionMat() creates the correct order 2 block companion matrix.",
      {
        testthat::skip_on_cran()
        beta <- list(
          matrix(
            data = c(
              0.40, 0.10,
              0.20, 0.30
            ),
            nrow = 2,
            ncol = 2
          ),
          matrix(
            data = c(
              0.10, 0.00,
              0.05, 0.10
            ),
            nrow = 2,
            ncol = 2
          )
        )

        expected <- matrix(
          data = c(
            0.40, 0.10, 1.00, 0.00,
            0.20, 0.30, 0.00, 1.00,
            0.10, 0.00, 0.00, 0.00,
            0.05, 0.10, 0.00, 0.00
          ),
          nrow = 4,
          ncol = 4
        )

        out <- CompanionMat(x = beta)

        testthat::expect_equal(
          object = out,
          expected = expected
        )
      }
    )

    testthat::test_that(
      "CompanionMat() creates the correct order 3 block companion matrix.",
      {
        testthat::skip_on_cran()
        beta <- list(
          matrix(
            data = c(
              0.10, 0.20,
              0.30, 0.40
            ),
            nrow = 2,
            ncol = 2
          ),
          matrix(
            data = c(
              0.50, 0.60,
              0.70, 0.80
            ),
            nrow = 2,
            ncol = 2
          ),
          matrix(
            data = c(
              0.90, 1.00,
              1.10, 1.20
            ),
            nrow = 2,
            ncol = 2
          )
        )

        expected <- matrix(
          data = 0.0,
          nrow = 6,
          ncol = 6
        )

        expected[1:2, 1:2] <- beta[[1]]
        expected[1:2, 3:4] <- beta[[2]]
        expected[1:2, 5:6] <- beta[[3]]
        expected[3:4, 1:2] <- diag(2)
        expected[5:6, 3:4] <- diag(2)

        out <- CompanionMat(x = beta)

        testthat::expect_equal(
          object = out,
          expected = expected
        )
      }
    )

    testthat::test_that(
      "CompanionMat() handles integer matrices by converting them to double.",
      {
        testthat::skip_on_cran()
        beta <- list(
          matrix(
            data = c(
              1L, 2L,
              3L, 4L
            ),
            nrow = 2,
            ncol = 2
          ),
          matrix(
            data = c(
              5L, 6L,
              7L, 8L
            ),
            nrow = 2,
            ncol = 2
          )
        )

        expected <- matrix(
          data = 0.0,
          nrow = 4,
          ncol = 4
        )

        expected[1:2, 1:2] <- beta[[1]]
        expected[1:2, 3:4] <- beta[[2]]
        expected[3:4, 1:2] <- diag(2)

        out <- CompanionMat(x = beta)

        testthat::expect_equal(
          object = out,
          expected = expected
        )
        testthat::expect_type(
          object = out,
          type = "double"
        )
      }
    )
  },
  text = "test-simStateSpace-companion-mat",
  tol = 0.01
)
