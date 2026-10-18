## ---- test-simStateSpace-check-simulation
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        p <- 2
        k <- 3
        j <- 2
        n <- 2
        time <- 4
        x <- lapply(
          X = seq_len(n),
          FUN = function(i) {
            matrix(data = 0.1 * seq_len(j * time), nrow = j)
          }
        )
        args <- list(
          n = n,
          time = time,
          mu0 = rep(x = 0, times = p),
          sigma0_l = diag(p),
          alpha = rep(x = 0, times = p),
          beta = diag(0.5, p),
          psi_l = diag(p),
          nu = rep(x = 0, times = k),
          lambda = matrix(data = 1, nrow = k, ncol = p),
          theta_l = diag(k),
          type = 2,
          x = x,
          gamma = matrix(data = 0.1, nrow = p, ncol = j),
          kappa = matrix(data = 0.1, nrow = k, ncol = j)
        )
        Run <- function(fun = SimSSMFixed, ...) {
          a <- args
          val <- list(...)
          a[names(val)] <- val
          do.call(what = fun, args = a)
        }
        Err <- function(regexp, ...) {
          testthat::expect_error(
            object = Run(...),
            regexp = regexp,
            fixed = TRUE
          )
        }
        # valid inputs pass, with vectors also as one-column matrices
        testthat::expect_s3_class(
          object = Run(alpha = matrix(data = 0, nrow = p)),
          class = "simstatespace"
        )
        Err("`n` must be a positive integer.", n = 1.5)
        Err("`time` must be a positive integer.", time = 0)
        Err("`delta_t` must be a positive number.", delta_t = -1)
        Err("`type` must be one of the types", type = 3)
        Err("`x` and `gamma` are required", gamma = NULL, type = 1)
        Err("`kappa` is required", kappa = NULL)
        Err(
          "`beta` must be a numeric 2 by 2 matrix.",
          beta = diag(0.5, p + 1)
        )
        Err("`alpha` must be a numeric vector of length 2.", alpha = 0)
        Err(
          "`theta_l` must be a numeric 3 by 3 matrix.",
          theta_l = diag(c(1, 1, NA))
        )
        Err(
          "`kappa` must be a numeric 3 by 2 matrix.",
          kappa = matrix(data = 0, nrow = k, ncol = j + 1)
        )
        Err(
          "`x` must be a list of numeric matrices with 2 rows",
          x = lapply(X = x, FUN = function(i) i[, -1])
        )
        # VAR: types 0 and 1 only
        testthat::expect_error(
          object = SimSSMVARFixed(
            n = n,
            time = time,
            mu0 = args$mu0,
            sigma0_l = args$sigma0_l,
            alpha = args$alpha,
            beta = args$beta,
            psi_l = args$psi_l,
            type = 2
          ),
          regexp = "`type` must be one of the types",
          fixed = TRUE
        )
        # individual-varying: every element is checked
        testthat::expect_error(
          object = SimSSMIVary(
            n = n,
            time = time,
            mu0 = list(args$mu0),
            sigma0_l = list(args$sigma0_l),
            alpha = list(args$alpha),
            beta = list(args$beta, diag(p + 1)),
            psi_l = list(args$psi_l),
            nu = list(args$nu),
            lambda = list(args$lambda),
            theta_l = list(args$theta_l)
          ),
          regexp = "Every element of `beta` must be a numeric 2 by 2 matrix.",
          fixed = TRUE
        )
        testthat::expect_error(
          object = SimSSMIVary(
            n = n,
            time = time,
            mu0 = args$mu0,
            sigma0_l = list(args$sigma0_l),
            alpha = list(args$alpha),
            beta = list(args$beta),
            psi_l = list(args$psi_l),
            nu = list(args$nu),
            lambda = list(args$lambda),
            theta_l = list(args$theta_l)
          ),
          regexp = "`mu0` must be a list with at least one element.",
          fixed = TRUE
        )
        # continuous-time parameters of the individual-varying linear SDE
        testthat::expect_error(
          object = SimSSMLinSDEIVary(
            n = n,
            time = time,
            mu0 = list(args$mu0),
            sigma0_l = list(args$sigma0_l),
            iota = list(c(0, 0, 0)),
            phi = list(-diag(p)),
            sigma_l = list(diag(p)),
            nu = list(args$nu),
            lambda = list(args$lambda),
            theta_l = list(args$theta_l)
          ),
          regexp = "Every element of `iota` must be a numeric vector of length",
          fixed = TRUE
        )
        testthat::expect_error(
          object = SimSSMOUIVary(
            n = n,
            time = time,
            mu0 = list(args$mu0),
            sigma0_l = list(args$sigma0_l),
            mu = list(args$mu0),
            phi = list(-diag(p)),
            sigma_l = list(diag(p + 1)),
            nu = list(args$nu),
            lambda = list(args$lambda),
            theta_l = list(args$theta_l)
          ),
          regexp = "Every element of `sigma_l` must be a numeric 2 by 2",
          fixed = TRUE
        )
        testthat::expect_error(
          object = SimSSMOUFixed(
            n = n,
            time = time,
            mu0 = args$mu0,
            sigma0_l = args$sigma0_l,
            mu = c(0, 0, 0),
            phi = -diag(p),
            sigma_l = diag(p),
            nu = args$nu,
            lambda = args$lambda,
            theta_l = args$theta_l
          ),
          regexp = "`mu` must be a numeric vector of length 2.",
          fixed = TRUE
        )
      }
    )
  },
  text = "test-simStateSpace-check-simulation"
)
