## ---- test-simStateSpace-wide-unbalanced
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        # Sim* output is balanced; drop one occasion of the second person
        # so the wide format falls back to stats::reshape()
        set.seed(42)
        sim <- simStateSpace::SimSSMOUFixed(
          n = 3,
          time = 4,
          delta_t = 0.5,
          mu0 = c(0, 0),
          sigma0_l = diag(x = 2),
          mu = c(0, 0),
          phi = -diag(x = 2),
          sigma_l = diag(x = 2),
          nu = c(0, 0),
          lambda = diag(x = 2),
          theta_l = diag(x = 2)
        )
        sim$data[[2]] <- lapply(
          X = sim$data[[2]],
          FUN = function(x) {
            x[-2, , drop = FALSE]
          }
        )
        long <- as.matrix(sim, eta = TRUE)
        wide <- as.matrix(
          stats::reshape(
            data = as.data.frame(long),
            timevar = "time",
            idvar = "id",
            direction = "wide",
            sep = "_"
          )
        )
        rownames(wide) <- NULL
        object <- as.matrix(sim, eta = TRUE, long = FALSE)
        attributes(object)[c("n", "m", "k", "p", "j")] <- NULL
        testthat::expect_identical(
          object = object,
          expected = wide
        )
        testthat::expect_true(
          is.na(object[2, "y1_0.5"])
        )
      }
    )
  },
  text = "test-simStateSpace-wide-unbalanced"
)
