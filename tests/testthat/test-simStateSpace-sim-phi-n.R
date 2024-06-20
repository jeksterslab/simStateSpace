## ---- test-simStateSpace-sim-phi-n
lapply(
  X = 1,
  FUN = function(i,
                 text,
                 n) {
    message(text)
    phi <- matrix(
      data = c(
        -0.357, 0.771, -0.450,
        0.0, -0.511, 0.729,
        0, 0, -0.693
      ),
      nrow = 3
    )
    vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
    testthat::test_that(
      text,
      {
        testthat::expect_true(
          all(
            phi - (1 / n) * Reduce(
              f = `+`,
              x = SimPhiN(
                n = n,
                phi = phi,
                vcov_phi_vec_l = vcov_phi_vec_l
              )
            ) < 0.001
          )
        )
      }
    )
  },
  text = "test-simStateSpace-sim-phi-n",
  n = 100000
)
