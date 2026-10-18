## ---- test-simStateSpace-transformations-reference
# Independent analytic fixtures complement the rounded chapter examples.
lapply(
  X = c(
    0.2,
    0.5,
    1.0,
    2.0
  ),
  FUN = function(delta_t,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      paste(text, delta_t),
      {
        testthat::skip_on_cran()
        # Scalar OU: evaluate closed-form integrals independently of C++.
        drift <- -0.4
        diffusion <- 1.7
        intercept <- 0.3
        beta <- exp(x = drift * delta_t)
        alpha <- intercept * expm1(x = drift * delta_t) / drift
        psi <- diffusion * expm1(x = 2 * drift * delta_t) / (2 * drift)
        out <- simStateSpace::LinSDE2SSM(
          iota = intercept,
          phi = matrix(data = drift),
          sigma_l = matrix(data = sqrt(x = diffusion)),
          delta_t = delta_t
        )
        testthat::expect_equal(
          object = c(out$alpha),
          expected = alpha,
          tolerance = tol
        )
        testthat::expect_equal(
          object = c(out$beta),
          expected = beta,
          tolerance = tol
        )
        testthat::expect_equal(
          object = c(tcrossprod(x = out$psi_l)),
          expected = psi,
          tolerance = tol
        )
        inverse <- simStateSpace::SSM2LinSDE(
          alpha = alpha,
          beta = matrix(data = beta),
          psi_l = matrix(data = sqrt(x = psi)),
          delta_t = delta_t
        )
        testthat::expect_equal(
          object = c(inverse$iota),
          expected = intercept,
          tolerance = tol
        )
        testthat::expect_equal(
          object = c(inverse$phi),
          expected = drift,
          tolerance = tol
        )
        testthat::expect_equal(
          object = c(tcrossprod(x = inverse$sigma_l)),
          expected = diffusion,
          tolerance = tol
        )
      }
    )
  },
  text = "test-simStateSpace-transformations-scalar-reference",
  tol = 1.0e-10
)

lapply(
  X = c(
    0.2,
    1.0,
    2.0
  ),
  FUN = function(delta_t,
                 text,
                 tol) {
    message(text)
    testthat::test_that(
      paste(text, delta_t),
      {
        testthat::skip_on_cran()
        # Brownian motion with drift: a previously singular solve.
        out <- simStateSpace::LinSDE2SSM(
          iota = c(
            0.3,
            -0.2
          ),
          phi = matrix(
            data = 0,
            nrow = 2,
            ncol = 2
          ),
          sigma_l = diag(
            x = c(
              0.5,
              0
            )
          ),
          delta_t = delta_t
        )
        testthat::expect_equal(
          object = c(out$alpha),
          expected = c(
            0.3,
            -0.2
          ) * delta_t,
          tolerance = tol
        )
        testthat::expect_equal(
          object = out$beta,
          expected = diag(x = 2),
          tolerance = tol
        )
        testthat::expect_equal(
          object = tcrossprod(x = out$psi_l),
          expected = diag(
            x = c(
              0.25,
              0
            )
          ) * delta_t,
          tolerance = tol
        )
        inverse <- simStateSpace::SSM2LinSDE(
          alpha = out$alpha,
          beta = out$beta,
          psi_l = out$psi_l,
          delta_t = delta_t
        )
        testthat::expect_equal(
          object = c(inverse$iota),
          expected = c(
            0.3,
            -0.2
          ),
          tolerance = tol
        )
        testthat::expect_equal(
          object = tcrossprod(x = inverse$sigma_l),
          expected = diag(
            x = c(
              0.25,
              0
            )
          ),
          tolerance = tol
        )
        # Integrated acceleration: state (position, velocity).
        h <- delta_t
        out <- simStateSpace::CTVAR2SSM(
          iota = 0.3,
          phi = list(
            matrix(data = 0),
            matrix(data = 0)
          ),
          sigma_l = matrix(data = 2),
          delta_t = h
        )
        expected_beta <- matrix(
          data = c(
            1,
            0,
            h,
            1
          ),
          nrow = 2
        )
        expected_psi <- 4 * matrix(
          data = c(
            h^3 / 3,
            h^2 / 2,
            h^2 / 2,
            h
          ),
          nrow = 2
        )
        testthat::expect_equal(
          object = out$beta,
          expected = expected_beta,
          tolerance = tol
        )
        testthat::expect_equal(
          object = c(out$alpha),
          expected = 0.3 * c(
            h^2 / 2,
            h
          ),
          tolerance = tol
        )
        testthat::expect_equal(
          object = tcrossprod(x = out$psi_l),
          expected = expected_psi,
          tolerance = tol
        )
        testthat::expect_equal(
          object = out$lambda,
          expected = matrix(
            data = c(
              1,
              0
            ),
            nrow = 1
          )
        )
        inverse <- simStateSpace::SSM2LinSDE(
          alpha = out$alpha,
          beta = out$beta,
          psi_l = out$psi_l,
          delta_t = h
        )
        testthat::expect_equal(
          object = c(inverse$iota),
          expected = c(
            0,
            0.3
          ),
          tolerance = tol
        )
        testthat::expect_equal(
          object = inverse$phi,
          expected = matrix(
            data = c(
              0,
              0,
              1,
              0
            ),
            nrow = 2
          ),
          tolerance = tol
        )
        testthat::expect_equal(
          object = tcrossprod(x = inverse$sigma_l),
          expected = diag(
            x = c(
              0,
              4
            )
          ),
          tolerance = tol
        )
      }
    )
  },
  text = "test-simStateSpace-transformations-singular-drift",
  tol = 1.0e-9
)

lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        # exp(A t) for this upper-triangular drift has an analytic expression.
        phi <- matrix(
          data = c(
            -0.3,
            0,
            0.7,
            -0.8
          ),
          nrow = 2
        )
        sigma <- matrix(
          data = c(
            1,
            0.2,
            0.2,
            0.6
          ),
          nrow = 2
        )
        transition <- function(t) {
          matrix(
            data = c(
              exp(x = -0.3 * t),
              0,
              1.4 * (exp(x = -0.3 * t) - exp(x = -0.8 * t)),
              exp(x = -0.8 * t)
            ),
            nrow = 2
          )
        }
        h <- 0.7
        expected <- matrix(
          data = 0,
          nrow = 2,
          ncol = 2
        )
        for (i in 1:2) {
          for (j in 1:2) {
            expected[i, j] <- integrate(
              f = Vectorize(
                FUN = function(t) {
                  e <- transition(t = t)
                  (e %*% sigma %*% t(x = e))[i, j]
                }
              ),
              lower = 0,
              upper = h,
              rel.tol = 1.0e-12
            )$value
          }
        }
        out <- simStateSpace::LinSDE2SSM(
          iota = c(
            0.2,
            -0.1
          ),
          phi = phi,
          sigma_l = t(x = chol(x = sigma)),
          delta_t = h
        )
        testthat::expect_equal(
          object = out$beta,
          expected = transition(t = h),
          tolerance = 1.0e-11
        )
        testthat::expect_equal(
          object = tcrossprod(x = out$psi_l),
          expected = expected,
          tolerance = 1.0e-11
        )
        # Semigroup identity tests independent interval composition.
        half <- simStateSpace::LinSDE2SSM(
          iota = c(
            0.2,
            -0.1
          ),
          phi = phi,
          sigma_l = t(x = chol(x = sigma)),
          delta_t = h / 2
        )
        testthat::expect_equal(
          object = out$beta,
          expected = half$beta %*% half$beta,
          tolerance = 1.0e-11
        )
        testthat::expect_equal(
          object = out$alpha,
          expected = half$alpha + half$beta %*% half$alpha,
          tolerance = 1.0e-11
        )
        testthat::expect_equal(
          object = tcrossprod(x = out$psi_l),
          expected = tcrossprod(x = half$psi_l) +
            half$beta %*% tcrossprod(x = half$psi_l) %*% t(
              x = half$beta
            ),
          tolerance = 1.0e-11
        )
      }
    )
  },
  text = "test-simStateSpace-transformations-noncommuting-reference"
)

lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        testthat::expect_error(
          object = simStateSpace::CTVAR2SSM(
            iota = 0,
            phi = list(matrix(data = 0)),
            sigma_l = matrix(data = 1),
            delta_t = 0
          ),
          regexp = "positive"
        )
        testthat::expect_error(
          object = simStateSpace::LinSDE2SSM(
            iota = numeric(),
            phi = matrix(
              data = numeric(),
              nrow = 0,
              ncol = 0
            ),
            sigma_l = matrix(
              data = numeric(),
              nrow = 0,
              ncol = 0
            )
          ),
          regexp = "positive dimensions"
        )
        testthat::expect_error(
          object = simStateSpace::SSM2LinSDE(
            alpha = 0,
            beta = matrix(data = 1),
            psi_l = matrix(data = 1),
            delta_t = 1,
            tol = NA
          ),
          regexp = "tol"
        )
        testthat::expect_error(
          object = simStateSpace::SSM2LinSDE(
            alpha = 0,
            beta = matrix(data = -1),
            psi_l = matrix(data = 1),
            delta_t = 1
          ),
          regexp = "imaginary"
        )
        # A valid DT covariance need not imply a PSD continuous-time diffusion.
        testthat::expect_error(
          object = simStateSpace::SSM2LinSDE(
            alpha = c(
              0,
              0
            ),
            beta = diag(
              x = c(
                0.9,
                0.1
              )
            ),
            psi_l = t(
              x = chol(
                x = matrix(
                  data = c(
                    1,
                    0.99,
                    0.99,
                    1
                  ),
                  nrow = 2
                )
              )
            ),
            delta_t = 1
          ),
          regexp = "positive semidefinite"
        )
        out <- simStateSpace::LinSDE2SSM(
          iota = c(
            0,
            0
          ),
          phi = matrix(
            data = 0,
            nrow = 2,
            ncol = 2
          ),
          sigma_l = diag(
            x = c(
              1.0e-8,
              0
            )
          ),
          delta_t = 1
        )
        testthat::expect_equal(
          object = tcrossprod(x = out$psi_l) / 1.0e-16,
          expected = diag(
            x = c(
              1,
              0
            )
          ),
          tolerance = 1.0e-12
        )
        rank_one <- matrix(
          data = c(
            1,
            2,
            0,
            0
          ),
          nrow = 2
        )
        out <- simStateSpace::LinSDE2SSM(
          iota = c(
            0,
            0
          ),
          phi = matrix(
            data = 0,
            nrow = 2,
            ncol = 2
          ),
          sigma_l = rank_one,
          delta_t = 1
        )
        testthat::expect_equal(
          object = tcrossprod(x = out$psi_l),
          expected = tcrossprod(x = rank_one),
          tolerance = 1.0e-12
        )
      }
    )
  },
  text = "test-simStateSpace-transformations-invalid-transforms"
)

lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        phi <- diag(
          x = c(
            -1000,
            -0.5
          )
        )
        sigma <- matrix(
          data = c(
            1,
            0.2,
            0.2,
            0.7
          ),
          nrow = 2
        )
        for (h in c(
          1,
          10
        )) {
          out <- simStateSpace::LinSDE2SSM(
            iota = c(
              0.3,
              -0.2
            ),
            phi = phi,
            sigma_l = t(x = chol(x = sigma)),
            delta_t = h
          )
          sums <- outer(
            X = diag(x = phi),
            Y = diag(x = phi),
            FUN = "+"
          )
          expected <- sigma * expm1(x = sums * h) / sums
          testthat::expect_equal(
            object = tcrossprod(x = out$psi_l),
            expected = expected,
            tolerance = 1.0e-10
          )
          testthat::expect_equal(
            object = c(out$alpha),
            expected = c(
              0.3,
              -0.2
            ) * expm1(x = diag(x = phi) * h) / diag(x = phi),
            tolerance = 1.0e-10
          )
          testthat::expect_equal(
            object = out$beta,
            expected = diag(x = exp(x = diag(x = phi) * h)),
            tolerance = 1.0e-10
          )
        }
      }
    )
  },
  text = "test-simStateSpace-transformations-strongly-decaying-drift"
)

lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        # Jordan block: exp(a I + N) = exp(a) (I + N) when N^2 = 0.
        a <- -0.4
        drift <- matrix(
          data = c(
            a,
            0,
            0.7,
            a
          ),
          nrow = 2
        )
        expected <- exp(x = a) * matrix(
          data = c(
            1,
            0,
            0.7,
            1
          ),
          nrow = 2
        )
        testthat::expect_equal(
          object = simStateSpace::ExpMat(x = drift),
          expected = expected,
          tolerance = 1.0e-12
        )
        testthat::expect_equal(
          object = simStateSpace::LogMat(x = expected),
          expected = drift,
          tolerance = 1.0e-10
        )
        testthat::expect_error(
          object = simStateSpace::LogMat(x = matrix(data = -1)),
          regexp = "imaginary"
        )
        testthat::expect_error(
          object = simStateSpace::LogMat(
            x = diag(x = 2),
            tol = Inf
          ),
          regexp = "finite"
        )
        testthat::expect_error(
          object = simStateSpace::ExpMat(
            x = matrix(
              data = 1,
              nrow = 2,
              ncol = 1
            )
          ),
          regexp = "square"
        )
      }
    )
  },
  text = "test-simStateSpace-transformations-matrix-functions"
)
