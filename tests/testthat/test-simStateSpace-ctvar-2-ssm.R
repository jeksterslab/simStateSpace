## ---- test-simStateSpace-ctvar-2-ssm
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      paste(text, "column matrix iota"),
      {
        testthat::skip_on_cran()
        phi <- list(
          matrix(data = c(-0.5, 0.1, 0, -0.4), nrow = 2),
          matrix(data = c(-0.3, 0, 0.05, -0.2), nrow = 2)
        )
        sigma_l <- diag(c(0.3, 0.2))
        testthat::expect_identical(
          object = simStateSpace::CTVAR2SSM(
            iota = matrix(data = c(0.1, -0.2), ncol = 1),
            phi = phi,
            sigma_l = sigma_l,
            delta_t = 0.5
          ),
          expected = simStateSpace::CTVAR2SSM(
            iota = c(0.1, -0.2),
            phi = phi,
            sigma_l = sigma_l,
            delta_t = 0.5
          )
        )
      }
    )
  },
  text = "test-simStateSpace-ctvar-2-ssm"
)
