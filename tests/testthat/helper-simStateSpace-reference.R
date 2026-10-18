## ---- helper-simStateSpace-reference
# Plain R references for the simulators.
# Definitions only: nothing is computed when the file is loaded.

# Exact discretization of the linear SDE
# d eta = (iota + phi eta) dt + sigma_l dW
# by eigendecomposition of phi (phi must be diagonalizable).
.RefDiscretize <- function(iota,
                           phi,
                           sigma,
                           delta_t) {
  eig <- eigen(phi)
  d <- eig$values
  v <- eig$vectors
  v_inv <- solve(v)
  # integral of exp(z * s) over [0, delta_t], elementwise
  int <- function(z) {
    out <- (exp(z * delta_t) - 1) / z
    out[abs(z) < 1e-12] <- delta_t
    out
  }
  beta <- Re(v %*% (exp(d * delta_t) * v_inv))
  alpha <- Re(v %*% (int(d) * (v_inv %*% iota)))
  sigma_tilde <- v_inv %*% sigma %*% t(v_inv)
  psi <- Re(v %*% (sigma_tilde * int(outer(d, d, "+"))) %*% t(v))
  psi <- (psi + t(psi)) / 2
  if (all(psi == 0)) {
    psi_l <- psi
  } else {
    psi_l <- t(chol(psi))
  }
  list(
    alpha = as.vector(alpha),
    beta = beta,
    psi = psi,
    psi_l = psi_l
  )
}

# Standard normal draws in the order the simulators use them.
# `arma::randn` follows `set.seed()` but does not reproduce `stats::rnorm()`,
# so the draws come from the package itself, with every parameter set
# so that the output is the noise:
# `eta` holds the state noise and `y` the measurement noise.
# Per person: initial state, initial measurement,
# then state and measurement noise at each later time point.
# With `k = NULL` (latent-only models), only state noise is drawn.
.RefNoise <- function(n,
                      time,
                      p,
                      k = NULL) {
  if (is.null(k)) {
    sim <- SimSSMVARFixed(
      n = n,
      time = time,
      mu0 = rep(x = 0, times = p),
      sigma0_l = diag(p),
      alpha = rep(x = 0, times = p),
      beta = matrix(data = 0, nrow = p, ncol = p),
      psi_l = diag(p)
    )
  } else {
    sim <- SimSSMFixed(
      n = n,
      time = time,
      mu0 = rep(x = 0, times = p),
      sigma0_l = diag(p),
      alpha = rep(x = 0, times = p),
      beta = matrix(data = 0, nrow = p, ncol = p),
      psi_l = diag(p),
      nu = rep(x = 0, times = k),
      lambda = matrix(data = 0, nrow = k, ncol = p),
      theta_l = diag(k)
    )
  }
  lapply(
    X = sim$data,
    FUN = function(x) {
      list(
        eta = x$eta,
        y = x$y
      )
    }
  )
}

# Parameter of person `i`:
# a list holds one value per person (recycled), anything else is fixed.
.RefPick <- function(x,
                     i) {
  if (is.list(x)) {
    x[[(i - 1) %% length(x) + 1]]
  } else {
    x
  }
}

# The state space recursion on the noise from `.RefNoise()`.
# Leave `nu`, `lambda`, and `theta_l` as `NULL` for latent-only models
# (`y` equals `eta`).
# Type 1: pass `x` and `gamma`. Type 2: also pass `kappa`.
.RefSim <- function(noise,
                    delta_t,
                    mu0,
                    sigma0_l,
                    alpha,
                    beta,
                    psi_l,
                    nu = NULL,
                    lambda = NULL,
                    theta_l = NULL,
                    x = NULL,
                    gamma = NULL,
                    kappa = NULL) {
  lapply(
    X = seq_along(noise),
    FUN = function(i) {
      z_eta <- noise[[i]]$eta
      z_y <- noise[[i]]$y
      time <- nrow(z_eta)
      p <- ncol(z_eta)
      mu0_i <- .RefPick(x = mu0, i = i)
      sigma0_l_i <- .RefPick(x = sigma0_l, i = i)
      alpha_i <- .RefPick(x = alpha, i = i)
      beta_i <- .RefPick(x = beta, i = i)
      psi_l_i <- .RefPick(x = psi_l, i = i)
      latent <- is.null(nu)
      if (!latent) {
        nu_i <- .RefPick(x = nu, i = i)
        lambda_i <- .RefPick(x = lambda, i = i)
        theta_l_i <- .RefPick(x = theta_l, i = i)
        k <- length(nu_i)
      } else {
        k <- p
      }
      if (!is.null(x)) {
        x_i <- x[[i]]
        gamma_i <- .RefPick(x = gamma, i = i)
      }
      if (!is.null(kappa)) {
        kappa_i <- .RefPick(x = kappa, i = i)
      }
      eta <- matrix(data = 0, nrow = time, ncol = p)
      y <- matrix(data = 0, nrow = time, ncol = k)
      for (t in seq_len(time)) {
        if (t == 1) {
          eta_t <- mu0_i + sigma0_l_i %*% z_eta[t, ]
        } else {
          eta_t <- alpha_i + beta_i %*% eta[t - 1, ] + psi_l_i %*% z_eta[t, ]
        }
        if (!is.null(x)) {
          eta_t <- eta_t + gamma_i %*% x_i[, t]
        }
        if (latent) {
          y_t <- eta_t
        } else {
          y_t <- nu_i + lambda_i %*% eta_t + theta_l_i %*% z_y[t, ]
          if (!is.null(kappa)) {
            y_t <- y_t + kappa_i %*% x_i[, t]
          }
        }
        eta[t, ] <- eta_t
        y[t, ] <- y_t
      }
      out <- list(
        id = rep(x = i, times = time),
        time = (seq_len(time) - 1) * delta_t,
        y = y,
        eta = eta
      )
      if (!is.null(x)) {
        out$x <- t(x_i)
      }
      out
    }
  )
}

# Compare the data of a `simstatespace` object with `.RefSim()`.
.RefExpectData <- function(object,
                           ref,
                           tol = 1e-12) {
  data <- lapply(
    X = object$data,
    FUN = function(x) {
      x$id <- as.vector(x$id)
      x$time <- as.vector(x$time)
      x
    }
  )
  testthat::expect_equal(
    object = data,
    expected = ref,
    tolerance = tol
  )
}

# Check the model description and the function name.
.RefExpectModel <- function(object,
                            model,
                            type,
                            fixed,
                            fun) {
  testthat::expect_identical(
    object = object$model,
    expected = list(
      model = model,
      covariates = type > 0,
      fixed = fixed,
      vary_i = !fixed
    )
  )
  testthat::expect_identical(
    object = object$fun,
    expected = fun
  )
}

# Drop the attributes the methods leave behind.
.RefStrip <- function(x) {
  attributes(x) <- attributes(x)[c("dim", "dimnames")]
  x
}

# Check the methods of a `simstatespace` object against `.RefSim()`:
# long and wide layouts, column names, `burnin`, and `reset_time`.
.RefExpectMethods <- function(object,
                              ref,
                              burnin = 2,
                              tol = 1e-12) {
  m <- length(ref[[1]]$time)
  times <- ref[[1]]$time
  k <- ncol(ref[[1]]$y)
  p <- ncol(ref[[1]]$eta)
  y_names <- paste0("y", seq_len(k))
  eta_names <- paste0("eta", seq_len(p))
  x_names <- NULL
  if (!is.null(ref[[1]]$x)) {
    x_names <- paste0("x", seq_len(ncol(ref[[1]]$x)))
  }
  long_eta <- do.call(
    what = "rbind",
    args = lapply(
      X = ref,
      FUN = function(r) {
        cbind(r$id, r$time, r$y, r$eta, r$x)
      }
    )
  )
  colnames(long_eta) <- c("id", "time", y_names, eta_names, x_names)
  long_y <- long_eta[, c("id", "time", y_names, x_names), drop = FALSE]
  # long
  testthat::expect_equal(
    object = .RefStrip(as.matrix(object, eta = TRUE)),
    expected = long_eta,
    tolerance = tol
  )
  testthat::expect_equal(
    object = .RefStrip(as.matrix(object)),
    expected = long_y,
    tolerance = tol
  )
  testthat::expect_equal(
    object = as.data.frame(object, eta = TRUE),
    expected = as.data.frame(long_eta),
    tolerance = tol
  )
  testthat::expect_equal(
    object = as.data.frame(object),
    expected = as.data.frame(long_y),
    tolerance = tol
  )
  # burnin
  keep <- long_eta[, "time"] %in% times[-seq_len(burnin)]
  burn <- long_eta[keep, , drop = FALSE]
  testthat::expect_equal(
    object = .RefStrip(
      as.matrix(object, eta = TRUE, burnin = burnin, reset_time = FALSE)
    ),
    expected = burn,
    tolerance = tol
  )
  burn[, "time"] <- burn[, "time"] - times[burnin + 1]
  testthat::expect_equal(
    object = .RefStrip(
      as.matrix(object, eta = TRUE, burnin = burnin, reset_time = TRUE)
    ),
    expected = burn,
    tolerance = tol
  )
  testthat::expect_equal(
    object = min(as.matrix(object, burnin = burnin)[, "time"]),
    expected = 0
  )
  testthat::expect_error(
    as.matrix(object, burnin = m),
    regexp = "`burnin` should not be greater than the measurement occasions."
  )
  # wide
  for (eta in c(FALSE, TRUE)) {
    if (eta) {
      long <- long_eta
    } else {
      long <- long_y
    }
    vars <- setdiff(colnames(long), c("id", "time"))
    wide <- do.call(
      what = "rbind",
      args = lapply(
        X = seq_along(ref),
        FUN = function(i) {
          c(i, as.vector(t(long[long[, "id"] == i, vars, drop = FALSE])))
        }
      )
    )
    colnames(wide) <- c(
      "id",
      paste(
        rep(x = vars, times = m),
        rep(x = times, each = length(vars)),
        sep = "_"
      )
    )
    testthat::expect_equal(
      object = .RefStrip(as.matrix(object, eta = eta, long = FALSE)),
      expected = wide,
      tolerance = tol
    )
    testthat::expect_equal(
      object = as.data.frame(object, eta = eta, long = FALSE),
      expected = as.data.frame(wide),
      tolerance = tol
    )
  }
  # print and plot
  testthat::expect_output(
    print(object),
    regexp = object$fun
  )
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  testthat::expect_no_error(
    plot(object, id = 1:2, time = times[1:3])
  )
  testthat::expect_no_error(
    plot(object, eta = TRUE)
  )
  testthat::expect_no_error(
    plot(object, burnin = burnin, reset_time = FALSE)
  )
}

# Standard normal draws for the random generators:
# `n` draws of length `q`, made in calls of length `q`,
# as the generators make them (the stream depends on the call length).
.RefStdNormal <- function(n,
                          q) {
  do.call(
    what = "rbind",
    args = SimMVN(
      n = n,
      mu = rep(x = 0, times = q),
      sigma_l = diag(q)
    )
  )
}

# Draw `n` values, one row of `z` per attempt:
# `make(i, z_row)` builds the candidate for draw `i`,
# and `accept(candidate)` keeps it.
# The number of attempts is stored in `attr(, "attempts")`.
.RefDraw <- function(n,
                     z,
                     make,
                     accept = function(x) TRUE) {
  out <- vector(mode = "list", length = n)
  row <- 0
  for (i in seq_len(n)) {
    repeat {
      row <- row + 1
      stopifnot(row <= nrow(z))
      candidate <- make(i, z[row, ])
      if (accept(candidate)) {
        out[[i]] <- candidate
        break
      }
    }
  }
  attr(out, "attempts") <- row
  out
}

# Plain R versions of the acceptance tests and projections.
.RefStationary <- function(x,
                           margin = 1) {
  all(Mod(eigen(x, only.values = TRUE)$values) < margin)
}

.RefTestPhi <- function(x,
                        margin = 0,
                        auto_ubound = 0) {
  all(Re(eigen(x, only.values = TRUE)$values) < margin) &&
    all(diag(x) <= auto_ubound)
}

.RefInBounds <- function(x,
                         lbound = NULL,
                         ubound = NULL) {
  ok <- TRUE
  if (!is.null(lbound)) {
    ok <- ok && all(x >= lbound, na.rm = TRUE)
  }
  if (!is.null(ubound)) {
    ok <- ok && all(x <= ubound, na.rm = TRUE)
  }
  ok
}

.RefProjectToStability <- function(x,
                                   margin = 0.98,
                                   tol = 1e-12) {
  rho <- max(Mod(eigen(x, only.values = TRUE)$values))
  if (rho < margin) {
    x
  } else {
    (margin / (rho + tol)) * x
  }
}

.RefProjectToHurwitz <- function(x,
                                 margin = 1e-3) {
  abscissa <- max(Re(eigen(x, only.values = TRUE)$values))
  if (abscissa >= -margin) {
    x - (abscissa + margin) * diag(nrow(x))
  } else {
    x
  }
}

# Spectral abscissa at most `-margin` (the projection lands on `-margin`).
.RefHurwitz <- function(x,
                        margin = 0) {
  max(Re(eigen(x, only.values = TRUE)$values)) <= -margin + 1e-10
}
