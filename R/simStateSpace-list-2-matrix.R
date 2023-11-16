#' Simulation Output to Matrix
#'
#' This function converts the output of
#' [simStateSpace::SimSSM0()],
#' [simStateSpace::SimSSMOU()],
#' [simStateSpace::SimSSMVAR()],
#' [simStateSpace::SimSSM0Fixed()],
#' [simStateSpace::SimSSMOUFixed()], or
#' [simStateSpace::SimSSMVARFixed()]
#' to a matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x R object.
#'   Output of
#'   [simStateSpace::SimSSM0()],
#'   [simStateSpace::SimSSMOU()],
#'   [simStateSpace::SimSSMVAR()],
#'   [simStateSpace::SimSSM0Fixed()],
#'   [simStateSpace::SimSSMOUFixed()], or
#'   [simStateSpace::SimSSMVARFixed()].
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
#' I <- diag(k)
#' I_sqrt <- chol(I)
#' null_vec <- rep(x = 0, times = k)
#' n <- 5
#' mu0 <- null_vec
#' sigma0_sqrt <- I_sqrt
#' alpha <- null_vec
#' beta <- diag(x = 0.50, nrow = k)
#' psi_sqrt <- I_sqrt
#' nu <- null_vec
#' lambda <- I
#' theta_sqrt <- chol(diag(x = 0.50, nrow = k))
#' time <- 50
#' burn_in <- 0
#'
#' # generate data
#' ssm <- SimSSM0(
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
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
#' ssm <- SimSSM0Fixed(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0_sqrt = sigma0_sqrt,
#'   alpha = alpha,
#'   beta = beta,
#'   psi_sqrt = psi_sqrt,
#'   nu = nu,
#'   lambda = lambda,
#'   theta_sqrt = theta_sqrt,
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
