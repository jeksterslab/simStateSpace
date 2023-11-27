#' Simulation Output to Matrix
#'
#' This function converts the output of
#' [simStateSpace::SimSSM()],
#' [simStateSpace::SimSSMOU()],
#' [simStateSpace::SimSSMVAR()],
#' [simStateSpace::SimSSMFixed()],
#' [simStateSpace::SimSSMOUFixed()],
#' [simStateSpace::SimSSMVARFixed()],
#' [simStateSpace::SimSSMVary()],
#' [simStateSpace::SimSSMOUVary()], or
#' [simStateSpace::SimSSMVARVary()]
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
#'   [simStateSpace::SimSSMVary()],
#'   [simStateSpace::SimSSMOUVary()], or
#'   [simStateSpace::SimSSMVARVary()].
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
#' # prepare parameters
#' set.seed(42)
#' k <- p <- 3
#' iden <- diag(k)
#' iden_sqrt <- chol(iden)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- null_vec
#' sigma0_sqrt <- iden_sqrt
#' alpha <- null_vec
#' beta <- diag(x = 0.50, nrow = k)
#' psi_sqrt <- iden_sqrt
#' nu <- null_vec
#' lambda <- iden
#' theta_sqrt <- chol(diag(x = 0.50, nrow = k))
#' time <- 50
#' burn_in <- 0
#'
#' # generate data
#' ssm <- SimSSM(
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   type = 0,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' # list to matrix
#' mat <- Sim2Matrix(ssm, long = TRUE)
#' str(mat)
#' head(mat)
#' mat <- Sim2Matrix(ssm, long = FALSE)
#' str(mat)
#' head(mat)
#'
#' # generate data
#' ssm <- SimSSMFixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
#'   type = 0,
#'   time = time,
#'   burn_in = burn_in
#' )
#'
#' # list to matrix
#' mat <- Sim2Matrix(ssm, long = TRUE)
#' str(mat)
#' head(mat)
#' mat <- Sim2Matrix(ssm, long = FALSE)
#' str(mat)
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
  if (is.list(first$x)) {
    covariates <- TRUE
    j <- dim(first$x[[1]])[2]
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
