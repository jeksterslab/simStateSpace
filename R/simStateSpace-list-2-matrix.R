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
#' mat <- Sim2Matrix(ssm)
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
#' mat <- Sim2Matrix(ssm)
#' str(mat)
#' head(mat)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace misc
#' @export
Sim2Matrix <- function(x,
                       eta = FALSE) {
  if (length(x) == 3) {
    if (
      identical(
        names(x),
        c("y", "eta", "time")
      )
    ) {
      n1 <- TRUE
    } else {
      n1 <- FALSE
    }
  } else {
    n1 <- FALSE
  }
  foo <- function(x,
                  eta,
                  n1) {
    k <- dim(x$y)[2]
    colnames(x$y) <- paste0("y", seq_len(k))
    colnames(x$time) <- "time"
    if (!n1) {
      colnames(x$id) <- "id"
    }
    if (!eta) {
      x$eta <- NULL
    } else {
      p <- dim(x$eta)[2]
      colnames(x$eta) <- paste0("eta", seq_len(p))
    }
    return(
      do.call(
        what = "cbind",
        args = x
      )
    )
  }
  if (n1) {
    return(
      foo(x, eta = eta, n1 = n1)
    )
  } else {
    return(
      do.call(
        what = "rbind",
        args = lapply(
          X = x,
          FUN = foo,
          eta = eta,
          n1 = n1
        )
      )
    )
  }
}
