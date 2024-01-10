#' Simulation Output to Matrix
#'
#' This function converts the output of
#' [simStateSpace::SimSSM()],
#' [simStateSpace::SimSSMOU()],
#' [simStateSpace::SimSSMVAR()],
#' [simStateSpace::SimSSMFixed()],
#' [simStateSpace::SimSSMOUFixed()],
#' [simStateSpace::SimSSMVARFixed()],
#' [simStateSpace::SimSSMIVary()],
#' [simStateSpace::SimSSMOUIVary()], or
#' [simStateSpace::SimSSMVARIVary()]
#' to a matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x R object.
#'   Output of
#'   [simStateSpace::SimSSM()],
#'   [simStateSpace::SimSSMOU()],
#'   [simStateSpace::SimSSMVAR()],
#'   [simStateSpace::SimSSMFixed()],
#'   [simStateSpace::SimSSMOUFixed()],
#'   [simStateSpace::SimSSMVARFixed()],
#'   [simStateSpace::SimSSMIVary()],
#'   [simStateSpace::SimSSMOUIVary()], or
#'   [simStateSpace::SimSSMVARIVary()].
#' @param eta Logical.
#'   If `eta = TRUE`, include `eta`.
#'   If `eta = FALSE`, exclude `eta`.
#' @param long Logical.
#'   If `long = TRUE`, use long format.
#'   If `long = FALSE`, use wide format.
#'
#' @return Returns a matrix of simulated data.
#'
#' @examples
#' # SimSSM Function -------------------------------------------------
#' # prepare parameters
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' null_vec <- rep(x = 0, times = k)
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
#' x <- matrix(
#'   data = rnorm(n = k * (time + burn_in)),
#'   ncol = k
#' )
#'
#' # Type 0
#' ssm <- SimSSM(
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
#' mat <- Sim2Matrix(ssm, eta = TRUE)
#' head(mat)
#' mat <- Sim2Matrix(ssm, eta = FALSE)
#' head(mat)
#' mat <- Sim2Matrix(ssm, eta = TRUE, long = FALSE)
#' head(mat)
#' mat <- Sim2Matrix(ssm, eta = FALSE, long = FALSE)
#' head(mat)
#'
#' # Type 1
#' ssm <- SimSSM(
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 1,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' mat <- Sim2Matrix(ssm, eta = TRUE)
#' head(mat)
#' mat <- Sim2Matrix(ssm, eta = FALSE)
#' head(mat)
#' mat <- Sim2Matrix(ssm, eta = TRUE, long = FALSE)
#' head(mat)
#' mat <- Sim2Matrix(ssm, eta = FALSE, long = FALSE)
#' head(mat)
#'
#' # Type 2
#' ssm <- SimSSM(
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   alpha = alpha,
#'   beta = beta,
#'   psi = psi,
#'   nu = nu,
#'   lambda = lambda,
#'   theta = theta,
#'   gamma_y = gamma_y,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 2,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' mat <- Sim2Matrix(ssm, eta = TRUE)
#' head(mat)
#' mat <- Sim2Matrix(ssm, eta = FALSE)
#' head(mat)
#' mat <- Sim2Matrix(ssm, eta = TRUE, long = FALSE)
#' head(mat)
#' mat <- Sim2Matrix(ssm, eta = FALSE, long = FALSE)
#' head(mat)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace misc
#' @export
Sim2Matrix <- function(x,
                       eta = FALSE,
                       long = TRUE) {
  if (length(x) == 5) {
    if (
      identical(
        names(x),
        c("y", "eta", "x", "time", "id")
      )
    ) {
      n1 <- TRUE
    } else {
      n1 <- FALSE
    }
  } else {
    n1 <- FALSE
  }
  if (n1) {
    first <- x
  } else {
    first <- x[[1]]
  }
  k <- dim(first$y)[2]
  y_names <- paste0("y", seq_len(k))
  if (length(y_names) == 1) {
    y_names <- "y"
  }
  p <- dim(first$eta)[2]
  eta_names <- paste0("eta", seq_len(p))
  if (length(eta_names) == 1) {
    eta_names <- "eta"
  }
  if (is.matrix(first$x)) {
    covariates <- TRUE
    j <- dim(first$x)[2]
    x_names <- paste0("x", seq_len(j))
  } else {
    covariates <- FALSE
  }
  foo <- function(x,
                  eta) {
    colnames(x$y) <- y_names
    colnames(x$time) <- "time"
    colnames(x$id) <- "id"
    if (covariates) {
      colnames(x$x) <- x_names
    } else {
      x$x <- NULL
    }
    if (!eta) {
      x$eta <- NULL
    } else {
      colnames(x$eta) <- eta_names
    }
    return(
      do.call(
        what = "cbind",
        args = x
      )
    )
  }
  if (n1) {
    out <- foo(
      x = x,
      eta = eta
    )
  } else {
    out <- do.call(
      what = "rbind",
      args = lapply(
        X = x,
        FUN = foo,
        eta = eta
      )
    )
  }
  if (long) {
    return(out)
  } else {
    out <- as.matrix(
      stats::reshape(
        data = as.data.frame(out),
        timevar = "time",
        idvar = "id",
        direction = "wide",
        sep = "_"
      )
    )
    rownames(out) <- seq_len(dim(out)[1])
    return(out)
  }
}
