.Long <- function(x, eta = FALSE) {
  # column names
  if (x$model$n1) {
    first <- x$data
  } else {
    first <- x$data[[1]]
  }
  obs <- first$y
  k <- dim(obs)[2]
  y_names <- paste0("y", seq_len(k))
  varnames <- c(
    "id",
    "time",
    y_names
  )
  lats <- first$eta
  p <- dim(lats)[2]
  eta_names <- paste0("eta", seq_len(p))
  varnames <- c(
    varnames,
    eta_names
  )
  if (x$model$covariates) {
    covs <- first$x
    j <- dim(covs)[2]
    x_names <- paste0("x", seq_len(j))
    varnames <- c(
      varnames,
      x_names
    )
  } else {
    j <- 0
  }
  if (x$model$n1) {
    out <- do.call(
      what = "cbind",
      args = x$data
    )
  } else {
    out <- lapply(
      X = x$data,
      FUN = function(x) {
        return(
          do.call(
            what = "cbind",
            args = x
          )
        )
      }
    )
    out <- do.call(
      what = "rbind",
      args = out
    )
  }
  colnames(out) <- varnames
  if (!eta) {
    varnames <- varnames[!(varnames %in% eta_names)]
    out <- out[, varnames, drop = FALSE]
  }
  attributes(out)$n <- length(
    unique(out[, "id"])
  )
  attributes(out)$k <- k
  attributes(out)$p <- p
  attributes(out)$j <- j
  return(out)
}

.Wide <- function(x, eta = FALSE) {
  long <- .Long(
    x = x,
    eta = eta
  )
  dims <- attributes(long)
  out <- as.matrix(
    stats::reshape(
      data = as.data.frame(
        long
      ),
      timevar = "time",
      idvar = "id",
      direction = "wide",
      sep = "_"
    )
  )
  rownames(out) <- seq_len(dim(out)[1])
  attributes(out)$n <- dims$n
  attributes(out)$k <- dims$k
  attributes(out)$p <- dims$p
  attributes(out)$j <- dims$j
  return(out)
}

#' Coerce an Object of Class `ssm` to a Data Frame
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Object of class `ssm`.
#' @param row.names `NULL` or character vector giving the row names
#'   for the data frame.
#'   Missing values are not allowed.
#' @param optional Logical.
#'   If `TRUE`, setting row names and converting column names is optional.
#' @param eta Logical.
#'   If `eta = TRUE`, include `eta`.
#'   If `eta = FALSE`, exclude `eta`.
#' @param long Logical.
#'   If `long = TRUE`, use long format.
#'   If `long = FALSE`, use wide format.
#' @param ... Additional arguments.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- null_vec
#' sigma0 <- iden
#' alpha <- null_vec
#' beta <- diag(x = 0.50, nrow = k)
#' psi <- iden
#' nu <- null_vec
#' lambda <- iden
#' theta <- diag(x = 0.50, nrow = k)
#' time <- 50
#' burn_in <- 0
#' gamma_y <- gamma_eta <- 0.10 * diag(k)
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = k * (time + burn_in)),
#'         ncol = k
#'       )
#'     )
#'   }
#' )
#'
#' # Type 0
#' ssm <- SimSSMFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   type = 0,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' as.data.frame(ssm)
#'
#' @keywords methods
#' @export
# nolint start: object_name_linter. Do not lint dot separators
as.data.frame.ssm <- function(x, row.names = NULL, # nolint
                              optional = FALSE,
                              eta = FALSE,
                              long = TRUE,
                              ...) {
  if (long) {
    out <- .Long(
      x = x,
      eta = eta
    )
  } else {
    out <- .Wide(
      x = x,
      eta = eta
    )
  }
  attributes(out)$n <- NULL
  attributes(out)$k <- NULL
  attributes(out)$p <- NULL
  attributes(out)$j <- NULL
  return(
    as.data.frame.matrix(
      x = out,
      row.names = row.names, # nolint
      optional = optional
    )
  )
}
# nolint end

#' Coerce an Object of Class `ssm` to a Matrix
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Object of class `ssm`.
#' @param eta Logical.
#'   If `eta = TRUE`, include `eta`.
#'   If `eta = FALSE`, exclude `eta`.
#' @param long Logical.
#'   If `long = TRUE`, use long format.
#'   If `long = FALSE`, use wide format.
#' @param ... Additional arguments.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- null_vec
#' sigma0 <- iden
#' alpha <- null_vec
#' beta <- diag(x = 0.50, nrow = k)
#' psi <- iden
#' nu <- null_vec
#' lambda <- iden
#' theta <- diag(x = 0.50, nrow = k)
#' time <- 50
#' burn_in <- 0
#' gamma_y <- gamma_eta <- 0.10 * diag(k)
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = k * (time + burn_in)),
#'         ncol = k
#'       )
#'     )
#'   }
#' )
#'
#' # Type 0
#' ssm <- SimSSMFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   type = 0,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' as.matrix(ssm)
#'
#' @keywords methods
#' @export
# nolint start: object_name_linter. Do not lint dot separators
as.matrix.ssm <- function(x,
                          eta = FALSE,
                          long = TRUE,
                          ...) {
  if (long) {
    out <- .Long(
      x = x,
      eta = eta
    )
  } else {
    out <- .Wide(
      x = x,
      eta = eta
    )
  }
  attributes(out)$n <- NULL
  attributes(out)$k <- NULL
  attributes(out)$p <- NULL
  attributes(out)$j <- NULL
  return(
    out
  )
}
# nolint end

#' Plot Method for an Object of Class `ssm`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Object of class `ssm`.
#' @param id Numeric vector.
#'   Optional `id` numbers to plot.
#'   If `id = NULL`, plot all available data.
#' @param time Numeric vector.
#'   Optional `time` points to plot.
#'   If `time = NULL`, plot all available data.
#' @param eta Logical.
#'   If `eta = TRUE`, plot the latent variables.
#'   If `eta = FALSE`, plot the observed variables.
#' @param type Character indicating the type of plotting;
#'   actually any of the types as in [plot.default()].
#' @param ... Additional arguments.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- null_vec
#' sigma0 <- iden
#' alpha <- null_vec
#' beta <- diag(x = 0.50, nrow = k)
#' psi <- iden
#' nu <- null_vec
#' lambda <- iden
#' theta <- diag(x = 0.50, nrow = k)
#' time <- 50
#' burn_in <- 0
#' gamma_y <- gamma_eta <- 0.10 * diag(k)
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = k * (time + burn_in)),
#'         ncol = k
#'       )
#'     )
#'   }
#' )
#'
#' # Type 0
#' ssm <- SimSSMFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   type = 0,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' plot(ssm)
#'
#' @keywords methods
#' @export
plot.ssm <- function(x,
                     id = NULL,
                     time = NULL,
                     eta = FALSE,
                     type = "b",
                     ...) {
  data <- .Long(
    x = x,
    eta = eta
  )
  if (eta) {
    n <- attributes(data)$p
    y <- paste0("eta", seq_len(n))
  } else {
    n <- attributes(data)$k
    y <- paste0("y", seq_len(n))
  }
  if (is.null(id)) {
    ids <- unique(data[, "id"])
  } else {
    ids <- id
  }
  colfunc <- grDevices::colorRampPalette(
    c(
      "red",
      "yellow",
      "springgreen",
      "royalblue"
    )
  )
  color <- colfunc(length(ids))
  for (i in seq_along(y)) {
    plot(
      x = 0,
      y = 0,
      xlim = range(data[, "time"]),
      ylim = range(data[, y]),
      type = "n",
      xlab = "time",
      ylab = y[i],
      main = y[i]
    )
    for (j in seq_along(ids)) {
      subset_data <- subset(
        x = data,
        subset = data[, "id"] == ids[j]
      )
      graphics::lines(
        x = subset_data[, "time"],
        y = subset_data[, y[i]],
        type = type,
        col = color[j],
        ...
      )
    }
  }
}

#' Print Method for an Object of Class `ssm`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Prints simulated data in long format.
#'
#' @param x Object of Class `ssm`.
#' @param eta Logical.
#'   If `eta = TRUE`, include `eta`.
#'   If `eta = FALSE`, exclude `eta`.
#' @param ... additional arguments.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- null_vec
#' sigma0 <- iden
#' alpha <- null_vec
#' beta <- diag(x = 0.50, nrow = k)
#' psi <- iden
#' nu <- null_vec
#' lambda <- iden
#' theta <- diag(x = 0.50, nrow = k)
#' time <- 50
#' burn_in <- 0
#' gamma_y <- gamma_eta <- 0.10 * diag(k)
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = k * (time + burn_in)),
#'         ncol = k
#'       )
#'     )
#'   }
#' )
#'
#' # Type 0
#' ssm <- SimSSMFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   type = 0,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' print(ssm)
#'
#' @keywords methods
#' @export
print.ssm <- function(x,
                      eta = FALSE,
                      ...) {
  out <- .Long(
    x = x,
    eta = eta
  )
  attributes(out)$n <- NULL
  attributes(out)$k <- NULL
  attributes(out)$p <- NULL
  attributes(out)$j <- NULL
  base::print(
    out
  )
}
