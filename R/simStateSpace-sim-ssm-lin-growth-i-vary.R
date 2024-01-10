#' Simulate Data from a Linear Growth Curve Model
#' (Individual-Varying Parameters)
#'
#' This function simulates data
#' from a linear growth curve model
#' for `n > 1` individuals.
#' In this model,
#' the parameters can vary across individuals.
#'
#' @details Parameters can vary across individuals
#'   by providing a list of parameter values.
#'   If the length of any of the parameters
#'   (`mu0`,
#'   `sigma0`,
#'   `mu`,
#'   `theta`,
#'   `gamma_y`, or
#'   `gamma_eta`)
#'   is less the `n`,
#'   the function will cycle through the available values.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param mu0 A list of numeric vectors.
#'   Each element of the list
#'   is a vector of length two.
#'   The first element is the mean of the intercept,
#'   and the second element is the mean of the slope.
#' @param sigma0 A list of numeric matrices.
#'   Each element of the list
#'   is the covariance matrix
#'   of the intercept and the slope.
#' @param theta A list numeric values.
#'   Each element of the list
#'   is the common measurement error variance.
#'
#' @inheritParams SimSSMLinGrowth
#' @inherit SimSSMFixed return
#' @inherit SimSSM references
#'
#' @examples
#' # prepare parameters
#' # In this example, the mean vector of the intercept and slope vary.
#' # Specifically,
#' # there are two sets of values representing two latent classes.
#' set.seed(42)
#' n <- 10
#' mu0_1 <- c(0.615, 1.006) # lower starting point, higher growth
#' mu0_2 <- c(1.000, 0.500) # higher starting point, lower growth
#' mu0 <- list(mu0_1, mu0_2)
#' sigma0 <- list(
#'   matrix(
#'     data = c(
#'       1.932,
#'       0.618,
#'       0.618,
#'       0.587
#'     ),
#'     nrow = 2
#'   )
#' )
#' theta <- list(0.6)
#' time <- 10
#' gamma_y <- list(matrix(data = 0.10, nrow = 1, ncol = 2))
#' gamma_eta <- list(matrix(data = 0.10, nrow = 2, ncol = 2))
#' x <- lapply(
#'   X = seq_len(n),
#'   FUN = function(i) {
#'     return(
#'       matrix(
#'         data = rnorm(n = 2 * time),
#'         ncol = 2
#'       )
#'     )
#'   }
#' )
#'
#' # Type 0
#' ssm <- SimSSMLinGrowthIVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   theta = theta,
#'   type = 0,
#'   time = time
#' )
#'
#' str(ssm)
#'
#' # Type 1
#' ssm <- SimSSMLinGrowthIVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   theta = theta,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 1,
#'   time = time
#' )
#'
#' str(ssm)
#'
#' # Type 2
#' ssm <- SimSSMLinGrowthIVary(
#'   n = n,
#'   mu0 = mu0,
#'   sigma0 = sigma0,
#'   theta = theta,
#'   gamma_y = gamma_y,
#'   gamma_eta = gamma_eta,
#'   x = x,
#'   type = 2,
#'   time = time
#' )
#'
#' str(ssm)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace sim growth
#' @export
SimSSMLinGrowthIVary <- function(n,
                                 mu0,
                                 sigma0,
                                 theta,
                                 gamma_y = NULL,
                                 gamma_eta = NULL,
                                 x = NULL,
                                 type = 0,
                                 time) {
  foo <- function(x) {
    return(
      t(chol(x))
    )
  }
  sigma0_l <- lapply(
    X = sigma0,
    FUN = foo
  )
  theta_l <- lapply(
    X = theta,
    FUN = sqrt
  )
  switch(
    EXPR = as.character(type),
    "0" = {
      return(
        .SimSSM0LinGrowthIVary(
          n = n,
          mu0 = rep(x = mu0, length.out = n),
          sigma0_l = rep(x = sigma0_l, length.out = n),
          theta_l = rep(x = theta_l, length.out = n),
          time = time
        )
      )
    },
    "1" = {
      return(
        .SimSSM1LinGrowthIVary(
          n = n,
          mu0 = rep(x = mu0, length.out = n),
          sigma0_l = rep(x = sigma0_l, length.out = n),
          theta_l = rep(x = theta_l, length.out = n),
          gamma_eta = rep(x = gamma_eta, length.out = n),
          x = x,
          time = time
        )
      )
    },
    "2" = {
      return(
        .SimSSM2LinGrowthIVary(
          n = n,
          mu0 = rep(x = mu0, length.out = n),
          sigma0_l = rep(x = sigma0_l, length.out = n),
          theta_l = rep(x = theta_l, length.out = n),
          gamma_y = rep(x = gamma_y, length.out = n),
          gamma_eta = rep(x = gamma_eta, length.out = n),
          x = x,
          time = time
        )
      )
    },
    stop(
      "Invalid `type`."
    )
  )
}
