#' Coerce an Object of Class `simstatespace` to a Data Frame
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Object of class `simstatespace`.
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
#' @param burnin Positive integer.
#'   Initial data points to discard.
#'   Default is zero.
#' @param reset_time Logical.
#'   Reset the time index after burnin.
#' @param ... Additional arguments.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' ## dynamic structure
#' p <- 3
#' mu0 <- rep(x = 0, times = p)
#' sigma0 <- diag(p)
#' sigma0_l <- t(chol(sigma0))
#' alpha <- rep(x = 0, times = p)
#' beta <- 0.50 * diag(p)
#' psi <- diag(p)
#' psi_l <- t(chol(psi))
#' ## measurement model
#' k <- 3
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta <- 0.50 * diag(k)
#' theta_l <- t(chol(theta))
#' ## covariates
#' j <- 2
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     matrix(
#'       data = stats::rnorm(n = time * j),
#'       nrow = j,
#'       ncol = time
#'     )
#'   }
#' )
#' gamma <- diag(x = 0.10, nrow = p, ncol = j)
#' kappa <- diag(x = 0.10, nrow = k, ncol = j)
#'
#' # Type 0
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 0
#' )
#'
#' head(as.data.frame(ssm))
#' head(as.data.frame(ssm, long = FALSE))
#'
#' # Type 1
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 1,
#'   x = x,
#'   gamma = gamma
#' )
#'
#' head(as.data.frame(ssm))
#' head(as.data.frame(ssm, long = FALSE))
#'
#' # Type 2
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 2,
#'   x = x,
#'   gamma = gamma,
#'   kappa = kappa
#' )
#'
#' head(as.data.frame(ssm))
#' head(as.data.frame(ssm, long = FALSE))
#'
#' @keywords methods
#' @export
as.data.frame.simstatespace <- function(x,
                                        row.names = NULL, # nolint: object_name_linter
                                        optional = FALSE,
                                        eta = FALSE,
                                        long = TRUE,
                                        burnin = 0,
                                        reset_time = TRUE,
                                        ...) {
  if (long) {
    out <- .Long(
      x = x,
      eta = eta,
      burnin = burnin,
      reset_time = reset_time
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
  as.data.frame.matrix(
    x = out,
    row.names = row.names,
    optional = optional
  )
}

#' Coerce an Object of Class `simstatespace` to a Matrix
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Object of class `simstatespace`.
#' @param eta Logical.
#'   If `eta = TRUE`, include `eta`.
#'   If `eta = FALSE`, exclude `eta`.
#' @param long Logical.
#'   If `long = TRUE`, use long format.
#'   If `long = FALSE`, use wide format.
#' @param burnin Positive integer.
#'   Initial data points to discard.
#'   Default is zero.
#' @param reset_time Logical.
#'   Reset the time index after burnin.
#' @param ... Additional arguments.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' ## dynamic structure
#' p <- 3
#' mu0 <- rep(x = 0, times = p)
#' sigma0 <- diag(p)
#' sigma0_l <- t(chol(sigma0))
#' alpha <- rep(x = 0, times = p)
#' beta <- 0.50 * diag(p)
#' psi <- diag(p)
#' psi_l <- t(chol(psi))
#' ## measurement model
#' k <- 3
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta <- 0.50 * diag(k)
#' theta_l <- t(chol(theta))
#' ## covariates
#' j <- 2
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     matrix(
#'       data = stats::rnorm(n = time * j),
#'       nrow = j,
#'       ncol = time
#'     )
#'   }
#' )
#' gamma <- diag(x = 0.10, nrow = p, ncol = j)
#' kappa <- diag(x = 0.10, nrow = k, ncol = j)
#'
#' # Type 0
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 0
#' )
#'
#' head(as.matrix(ssm))
#' head(as.matrix(ssm, long = FALSE))
#'
#' # Type 1
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 1,
#'   x = x,
#'   gamma = gamma
#' )
#'
#' head(as.matrix(ssm))
#' head(as.matrix(ssm, long = FALSE))
#'
#' # Type 2
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 2,
#'   x = x,
#'   gamma = gamma,
#'   kappa = kappa
#' )
#'
#' head(as.matrix(ssm))
#' head(as.matrix(ssm, long = FALSE))
#'
#' @keywords methods
#' @export
as.matrix.simstatespace <- function(x,
                                    eta = FALSE,
                                    long = TRUE,
                                    burnin = 0,
                                    reset_time = TRUE,
                                    ...) {
  if (long) {
    out <- .Long(
      x = x,
      eta = eta,
      burnin = burnin,
      reset_time = reset_time
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
  out
}

#' Plot Method for an Object of Class `simstatespace`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x Object of class `simstatespace`.
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
#' @param burnin Positive integer.
#'   Initial data points to discard.
#'   Default is zero.
#' @param reset_time Logical.
#'   Reset the time index after burnin.
#' @param ... Additional arguments.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' ## dynamic structure
#' p <- 3
#' mu0 <- rep(x = 0, times = p)
#' sigma0 <- diag(p)
#' sigma0_l <- t(chol(sigma0))
#' alpha <- rep(x = 0, times = p)
#' beta <- 0.50 * diag(p)
#' psi <- diag(p)
#' psi_l <- t(chol(psi))
#' ## measurement model
#' k <- 3
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta <- 0.50 * diag(k)
#' theta_l <- t(chol(theta))
#' ## covariates
#' j <- 2
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     matrix(
#'       data = stats::rnorm(n = time * j),
#'       nrow = j,
#'       ncol = time
#'     )
#'   }
#' )
#' gamma <- diag(x = 0.10, nrow = p, ncol = j)
#' kappa <- diag(x = 0.10, nrow = k, ncol = j)
#'
#' # Type 0
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 0
#' )
#'
#' plot(ssm)
#' plot(ssm, id = 1:3, time = 0:9)
#'
#' # Type 1
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 1,
#'   x = x,
#'   gamma = gamma
#' )
#'
#' plot(ssm)
#' plot(ssm, id = 1:3, time = 0:9)
#'
#' # Type 2
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 2,
#'   x = x,
#'   gamma = gamma,
#'   kappa = kappa
#' )
#'
#' plot(ssm)
#' plot(ssm, id = 1:3, time = 0:9)
#'
#' @keywords methods
#' @export
plot.simstatespace <- function(x,
                               id = NULL,
                               time = NULL,
                               eta = FALSE,
                               type = "b",
                               burnin = 0,
                               reset_time = TRUE,
                               ...) {
  data <- .Long(
    x = x,
    eta = eta,
    burnin = burnin,
    reset_time = reset_time
  )
  if (eta) {
    n <- attributes(data)$p
    y <- paste0("eta", seq_len(n))
  } else {
    n <- attributes(data)$k
    y <- paste0("y", seq_len(n))
  }
  if (!is.null(id)) {
    data <- data[which(data[, "id"] %in% id), , drop = FALSE]
  }
  if (!is.null(time)) {
    data <- data[which(data[, "time"] %in% time), , drop = FALSE]
  }
  colfunc <- grDevices::colorRampPalette(
    c(
      "red",
      "yellow",
      "springgreen",
      "royalblue"
    )
  )
  ids <- unique(data[, "id"])
  color <- colfunc(length(ids))
  for (i in seq_along(y)) {
    graphics::plot.default(
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

#' Print Method for an Object of Class `simstatespace`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Prints simulated data in long format.
#'
#' @param x Object of Class `simstatespace`.
#' @param ... Additional arguments.
#'
#' @examples
#' # prepare parameters
#' set.seed(42)
#' ## number of individuals
#' n <- 5
#' ## time points
#' time <- 50
#' ## dynamic structure
#' p <- 3
#' mu0 <- rep(x = 0, times = p)
#' sigma0 <- diag(p)
#' sigma0_l <- t(chol(sigma0))
#' alpha <- rep(x = 0, times = p)
#' beta <- 0.50 * diag(p)
#' psi <- diag(p)
#' psi_l <- t(chol(psi))
#' ## measurement model
#' k <- 3
#' nu <- rep(x = 0, times = k)
#' lambda <- diag(k)
#' theta <- 0.50 * diag(k)
#' theta_l <- t(chol(theta))
#' ## covariates
#' j <- 2
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     matrix(
#'       data = stats::rnorm(n = time * j),
#'       nrow = j,
#'       ncol = time
#'     )
#'   }
#' )
#' gamma <- diag(x = 0.10, nrow = p, ncol = j)
#' kappa <- diag(x = 0.10, nrow = k, ncol = j)
#'
#' # Type 0
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 0
#' )
#'
#' print(ssm)
#'
#' # Type 1
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 1,
#'   x = x,
#'   gamma = gamma
#' )
#'
#' print(ssm)
#'
#' # Type 2
#' ssm <- SimSSMFixed(
#'   n = n,
#'   time = time,
#'   mu0 = mu0,
#'   sigma0_l = sigma0_l,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_l = psi_l,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_l = theta_l,
#'   type = 2,
#'   x = x,
#'   gamma = gamma,
#'   kappa = kappa
#' )
#'
#' print(ssm)
#'
#' @keywords methods
#' @export
print.simstatespace <- function(x,
                                ...) {
  cat("Call:\n")
  base::print(x$call)
  cat(
    paste0(
      "Use `as.data.frame` or `as.matrix` on the output of",
      " `",
      x$fun,
      "`",
      "\nto convert it to a data frame or a matrix.\n",
      "\n"
    )
  )
}
