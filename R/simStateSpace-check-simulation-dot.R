.CheckSimulation <- function(n, time, delta_t,
                             par,
                             type, types = 0:2,
                             x = NULL, gamma = NULL, kappa = NULL,
                             fixed = TRUE) {
  # Input checks shared by every `Sim*` function.
  # `par` is a named list of the model parameters other than `gamma` and
  # `kappa`: one set each with `fixed = TRUE`, lists of sets otherwise.
  # Vectors: `mu0`, `alpha`, `iota`, `mu`, `nu`; matrices: the rest.
  # Only structure is checked (numeric, finite, dimensions),
  # so estimated parameters always pass.
  Count <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x == round(x)
  }
  stopifnot(
    "`n` must be a positive integer." = Count(n),
    "`time` must be a positive integer." = Count(time),
    "`delta_t` must be a positive number." = is.numeric(delta_t) &&
      length(delta_t) == 1L && is.finite(delta_t) && delta_t > 0,
    "`type` must be one of the types this model supports." =
      length(type) == 1L && type %in% types
  )
  if (type > 0) {
    stopifnot(
      "`x` and `gamma` are required when `type` is 1 or 2." =
        !is.null(x) && !is.null(gamma)
    )
    par$gamma <- gamma
  }
  if (type > 1) {
    stopifnot(
      "`kappa` is required when `type` is 2." = !is.null(kappa)
    )
    par$kappa <- kappa
  }
  par <- par[!vapply(X = par, FUN = is.null, FUN.VALUE = logical(1))]
  if (fixed) {
    par <- lapply(X = par, FUN = list)
  }
  for (name in names(par)) {
    if (!is.list(par[[name]]) || length(par[[name]]) == 0L) {
      stop(
        "`", name, "` must be a list with at least one element.",
        call. = FALSE
      )
    }
  }
  p <- length(par$mu0[[1]])
  k <- if (is.null(par$nu)) p else length(par$nu[[1]])
  j <- if (is.null(par$gamma)) 0L else NCOL(par$gamma[[1]])
  dims <- list(
    mu0 = p, alpha = p, iota = p, mu = p, nu = k,
    sigma0_l = c(p, p), beta = c(p, p), psi_l = c(p, p),
    phi = c(p, p), sigma_l = c(p, p),
    lambda = c(k, p), theta_l = c(k, k),
    gamma = c(p, j), kappa = c(k, j)
  )
  for (name in names(par)) {
    d <- dims[[name]]
    ok <- vapply(
      X = par[[name]],
      FUN = function(i) {
        is.numeric(i) && all(is.finite(i)) && if (length(d) == 1L) {
          length(i) == d
        } else {
          is.matrix(i) && identical(x = dim(i), y = as.integer(d))
        }
      },
      FUN.VALUE = logical(1)
    )
    if (!all(ok)) {
      shape <- if (length(d) == 1L) {
        paste("a numeric vector of length", d)
      } else {
        paste0("a numeric ", d[1], " by ", d[2], " matrix")
      }
      stop(
        if (fixed) {
          paste0("`", name, "` must be ", shape, ".")
        } else {
          paste0("Every element of `", name, "` must be ", shape, ".")
        },
        call. = FALSE
      )
    }
  }
  if (type > 0) {
    ok <- is.list(x) && length(x) > 0L && all(
      vapply(
        X = x,
        FUN = function(i) {
          is.matrix(i) && is.numeric(i) && all(is.finite(i)) &&
            nrow(i) == j && ncol(i) == time
        },
        FUN.VALUE = logical(1)
      )
    )
    if (!ok) {
      stop(
        "`x` must be a list of numeric matrices with ", j,
        " rows (the columns of `gamma`) and `time` columns.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}
