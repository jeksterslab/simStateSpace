// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-var.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Data from the Vector Autoregressive Model
//' using a State Space Model Parameterization (n = 1)
//'
//' This function simulates data from the vector autoregressive model
//' using a state space model parameterization.
//' See details for more information.
//'
//' @details The measurement model is given by
//'   \deqn{
//'     \mathbf{y}_{t}
//'     =
//'     \boldsymbol{\eta}_{t} .
//'   }
//'
//'   The dynamic structure is given by
//'   \deqn{
//'     \boldsymbol{\eta}_{t}
//'     =
//'     \boldsymbol{\alpha}
//'     +
//'     \boldsymbol{\beta}
//'     \boldsymbol{\eta}_{t - 1}
//'     +
//'     \boldsymbol{\zeta}_{t}
//'     \quad
//'     \mathrm{with}
//'     \quad
//'     \boldsymbol{\zeta}_{t}
//'     \sim
//'     \mathcal{N}
//'     \left(
//'     \mathbf{0},
//'     \boldsymbol{\Psi}
//'     \right)
//'   }
//'   where \eqn{\boldsymbol{\eta}_{t}}, \eqn{\boldsymbol{\eta}_{t - 1}},
//'   and \eqn{\boldsymbol{\zeta}_{t}} are random variables
//'   and \eqn{\boldsymbol{\alpha}}, \eqn{\boldsymbol{\beta}},
//'   and \eqn{\boldsymbol{\Psi}} are model parameters.
//'   \eqn{\boldsymbol{\eta}_{t}} is a vector of latent variables
//'   at time \eqn{t}, \eqn{\boldsymbol{\eta}_{t - 1}}
//'   is a vector of latent variables at
//'   \eqn{t - 1},
//'   and \eqn{\boldsymbol{\zeta}_{t}} is a vector of dynamic noise
//'   at time \eqn{t} while \eqn{\boldsymbol{\alpha}}
//'   is a vector of intercepts,
//'   \eqn{\boldsymbol{\beta}} is a matrix of autoregression
//'   and cross regression coefficients,
//'   and \eqn{\boldsymbol{\Psi}} is the covariance matrix of
//'   \eqn{\boldsymbol{\zeta}_{t}}.
//'
//' @inheritParams SimSSM0
//' @inherit SimSSM0 return
//' @inherit SimSSM0 references
//'
//' @examples
//' # prepare parameters
//' set.seed(42)
//' k <- 3
//' I <- diag(k)
//' I_sqrt <- chol(I)
//' null_vec <- rep(x = 0, times = k)
//' mu0 <- null_vec
//' sigma0_sqrt <- I_sqrt
//' alpha <- null_vec
//' beta <- diag(x = 0.5, nrow = k)
//' psi_sqrt <- I_sqrt
//' time <- 50
//' burn_in <- 0
//'
//' # generate data
//' ssm <- SimSSMVAR(
//'   mu0 = mu0,
//'   sigma0_sqrt = sigma0_sqrt,
//'   alpha = alpha,
//'   beta = beta,
//'   psi_sqrt = psi_sqrt,
//'   time = time,
//'   burn_in = burn_in
//' )
//'
//' str(ssm)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim
//' @export
// [[Rcpp::export]]
Rcpp::List SimSSMVAR(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                     const arma::vec& alpha, const arma::mat& beta,
                     const arma::mat& psi_sqrt, const int time,
                     const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);

  // Step 3: Generate initial condition
  eta.col(0) = mu0 + sigma0_sqrt * arma::randn(num_latent_vars);

  // Step 4: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) =
        alpha + beta * eta.col(t - 1) + psi_sqrt * arma::randn(num_latent_vars);
  }

  // Step 5: If there is a burn-in period, remove it
  if (burn_in > 0) {
    eta = eta.cols(burn_in, total_time - 1);
  }

  // Step 6: Return the transposed data matrices in a list
  return Rcpp::List::create(Rcpp::Named("y") = eta.t(),
                            Rcpp::Named("eta") = eta.t(),
                            Rcpp::Named("time") = arma::regspace(1, time));
}
