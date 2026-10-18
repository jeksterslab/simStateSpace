## ---- test-simStateSpace-sim-ssm-lin-sde-edge
lapply(
  X = 1,
  FUN = function(i,
                 text) {
    message(text)
    testthat::test_that(
      text,
      {
        testthat::skip_on_cran()
        # one edge case per person:
        # 1. singular drift (eigenvalues 0 and -0.5) with a nonzero iota
        # 2. semidefinite diffusion (noise on the first state only)
        # 3. undamped oscillator (eigenvalues +/- i, lambda_1 + lambda_2 = 0)
        n <- 3
        time <- 6
        delta_t <- 0.5
        p <- 2
        k <- 2
        mu0 <- list(c(1, -1), c(0, 0.5), c(0.5, 0))
        sigma0_l <- list(diag(c(0.5, 0.2)), diag(c(0.3, 0.4)), diag(2))
        iota <- list(c(0.2, -0.1), c(0.1, 0.3), c(0, 0))
        phi <- list(
          matrix(data = c(0, 0.3, 0, -0.5), nrow = p),
          diag(c(-0.5, -1)),
          matrix(data = c(0, -1, 1, 0), nrow = p)
        )
        sigma_l <- list(
          t(chol(matrix(data = c(0.3, 0.1, 0.1, 0.2), nrow = p))),
          matrix(data = c(0.6, 0, 0, 0), nrow = p),
          diag(c(0.4, 0.3))
        )
        nu <- c(0.5, -0.5)
        lambda <- diag(p)
        theta_l <- diag(c(0.2, 0.3))
        disc <- lapply(
          X = c(1, 3),
          FUN = function(i) {
            .RefDiscretize(
              iota = iota[[i]],
              phi = phi[[i]],
              sigma = sigma_l[[i]] %*% t(sigma_l[[i]]),
              delta_t = delta_t
            )
          }
        )
        # case 2 in closed form: psi = diag(0.36 * (1 - exp(-delta_t)), 0)
        disc <- list(
          disc[[1]],
          list(
            alpha = (exp(-c(0.5, 1) * delta_t) - 1) / -c(0.5, 1) * iota[[2]],
            beta = diag(exp(-c(0.5, 1) * delta_t)),
            psi_l = diag(c(sqrt(0.36 * (1 - exp(-delta_t))), 0))
          ),
          disc[[2]]
        )
        Ref <- function(d = disc, i = NULL, persons = n) {
          pick <- function(name) {
            out <- lapply(X = d, FUN = "[[", name)
            if (is.null(i)) out else out[[i]]
          }
          noise <- .RefNoise(n = persons, time = time, p = p, k = k)
          .RefSim(
            noise = noise,
            delta_t = delta_t,
            mu0 = if (is.null(i)) mu0 else mu0[[i]],
            sigma0_l = if (is.null(i)) sigma0_l else sigma0_l[[i]],
            alpha = pick("alpha"),
            beta = pick("beta"),
            psi_l = pick("psi_l"),
            nu = nu,
            lambda = lambda,
            theta_l = theta_l
          )
        }
        # individual-varying
        set.seed(42)
        ssm <- SimSSMLinSDEIVary(
          n = n,
          time = time,
          delta_t = delta_t,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          iota = iota,
          phi = phi,
          sigma_l = sigma_l,
          nu = list(nu),
          lambda = list(lambda),
          theta_l = list(theta_l),
          type = 0
        )
        set.seed(42)
        .RefExpectData(object = ssm, ref = Ref())
        # OU parameterization of the same models: iota = -phi mu;
        # case 1 (singular drift) has no mu, so use iota = 0 there
        mu <- list(
          c(0, 0),
          -solve(phi[[2]], iota[[2]]),
          c(0.4, -0.2)
        )
        disc_ou <- disc
        disc_ou[[1]]$alpha <- c(0, 0)
        disc_ou[[3]]$alpha <- .RefDiscretize(
          iota = -phi[[3]] %*% mu[[3]],
          phi = phi[[3]],
          sigma = sigma_l[[3]] %*% t(sigma_l[[3]]),
          delta_t = delta_t
        )$alpha
        set.seed(42)
        ssm <- SimSSMOUIVary(
          n = n,
          time = time,
          delta_t = delta_t,
          mu0 = mu0,
          sigma0_l = sigma0_l,
          mu = mu,
          phi = phi,
          sigma_l = sigma_l,
          nu = list(nu),
          lambda = list(lambda),
          theta_l = list(theta_l),
          type = 0
        )
        set.seed(42)
        .RefExpectData(object = ssm, ref = Ref(d = disc_ou))
        # fixed
        for (i in seq_len(n)) {
          set.seed(i)
          ssm <- SimSSMLinSDEFixed(
            n = 2,
            time = time,
            delta_t = delta_t,
            mu0 = mu0[[i]],
            sigma0_l = sigma0_l[[i]],
            iota = iota[[i]],
            phi = phi[[i]],
            sigma_l = sigma_l[[i]],
            nu = nu,
            lambda = lambda,
            theta_l = theta_l,
            type = 0
          )
          set.seed(i)
          .RefExpectData(object = ssm, ref = Ref(i = i, persons = 2))
        }
      }
    )
  },
  text = "test-simStateSpace-sim-ssm-lin-sde-edge"
)
