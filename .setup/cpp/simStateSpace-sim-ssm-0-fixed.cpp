// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-fixed.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Data using a State Space Model Parameterization
//' for n > 1 Individuals (Fixed Parameters)
//'
//' This function simulates data
//' using a state space model parameterization
//' for `n > 1` individuals.
//' In this model,
//' the parameters are invariant across individuals.
//'
//' @details The measurement model is given by
//'   \deqn{
//'     \mathbf{y}_{i, t}
//'     =
//'     \boldsymbol{\nu}
//'     +
//'     \boldsymbol{\Lambda}
//'     \boldsymbol{\eta}_{i, t}
//'     +
//'     \boldsymbol{\varepsilon}_{i, t}
//'     \quad
//'     \mathrm{with}
//'     \quad
//'     \boldsymbol{\varepsilon}_{i, t}
//'     \sim
//'     \mathcal{N}
//'     \left(
//'     \mathbf{0},
//'     \boldsymbol{\Theta}
//'     \right)
//'   }
//'   where \eqn{\mathbf{y}_{i, t}}, \eqn{\boldsymbol{\eta}_{i, t}},
//'   and \eqn{\boldsymbol{\varepsilon}_{i, t}}
//'   are random variables and \eqn{\boldsymbol{\nu}},
//'   \eqn{\boldsymbol{\Lambda}},
//'   and \eqn{\boldsymbol{\Theta}} are model parameters.
//'   \eqn{\mathbf{y}_{i, t}} is a vector of observed random variables
//'   at time \eqn{t} and individual \eqn{i},
//'   \eqn{\boldsymbol{\eta}_{i, t}} is a vector of latent random variables
//'   at time \eqn{t} and individual \eqn{i},
//'   and \eqn{\boldsymbol{\varepsilon}_{i, t}}
//'   is a vector of random measurement errors
//'   at time \eqn{t} and individual \eqn{i},
//'   while \eqn{\boldsymbol{\nu}} is a vector of intercept,
//'   \eqn{\boldsymbol{\Lambda}} is a matrix of factor loadings,
//'   and \eqn{\boldsymbol{\Theta}} is the covariance matrix of
//'   \eqn{\boldsymbol{\varepsilon}}.
//'
//'   The dynamic structure is given by
//'   \deqn{
//'     \boldsymbol{\eta}_{i, t}
//'     =
//'     \boldsymbol{\alpha}
//'     +
//'     \boldsymbol{\beta}
//'     \boldsymbol{\eta}_{i, t - 1}
//'     +
//'     \boldsymbol{\zeta}_{i, t}
//'     \quad
//'     \mathrm{with}
//'     \quad
//'     \boldsymbol{\zeta}_{i, t}
//'     \sim
//'     \mathcal{N}
//'     \left(
//'     \mathbf{0},
//'     \boldsymbol{\Psi}
//'     \right)
//'   }
//'   where \eqn{\boldsymbol{\eta}_{i, t}}, \eqn{\boldsymbol{\eta}_{i, t - 1}},
//'   and \eqn{\boldsymbol{\zeta}_{i, t}} are random variables
//'   and \eqn{\boldsymbol{\alpha}}, \eqn{\boldsymbol{\beta}},
//'   and \eqn{\boldsymbol{\Psi}} are model parameters.
//'   \eqn{\boldsymbol{\eta}_{i, t}} is a vector of latent variables
//'   at time \eqn{t} and individual \eqn{i},
//'   \eqn{\boldsymbol{\eta}_{i, t - 1}}
//'   is a vector of latent variables at
//'   time \eqn{t - 1} and individual \eqn{i},
//'   and \eqn{\boldsymbol{\zeta}_{i, t}} is a vector of dynamic noise
//'   at time \eqn{t} and individual \eqn{i} while \eqn{\boldsymbol{\alpha}}
//'   is a vector of intercepts,
//'   \eqn{\boldsymbol{\beta}} is a matrix of autoregression
//'   and cross regression coefficients,
//'   and \eqn{\boldsymbol{\Psi}} is the covariance matrix of
//'   \eqn{\boldsymbol{\zeta}_{i, t}}.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param n Positive integer.
//'   Number of individuals.
//' @inheritParams SimSSM0
//' @inherit SimSSM0 references
//'
//' @return Returns a list of length `n`.
//'   Each element is a list with the following elements:
//'   - `y`: A `t` by `k` matrix of values for the manifest variables.
//'   - `eta`: A `t` by `p` matrix of values for the latent variables.
//'   - `time`: A vector of discrete time points from 1 to `t`.
//'   - `id`: A vector of ID numbers of length `t`.
//'   - `n`: Number of individuals.
//'
//' @examples
//' # prepare parameters
//' set.seed(42)
//' k <- p <- 3
//' I <- diag(k)
//' I_sqrt <- chol(I)
//' null_vec <- rep(x = 0, times = k)
//' n <- 5
//' mu0 <- null_vec
//' sigma0_sqrt <- I_sqrt
//' alpha <- null_vec
//' beta <- diag(x = 0.50, nrow = k)
//' psi_sqrt <- I_sqrt
//' nu <- null_vec
//' lambda <- I
//' theta_sqrt <- chol(diag(x = 0.50, nrow = k))
//' time <- 50
//' burn_in <- 0
//'
//' # generate data
//' ssm <- SimSSM0Fixed(
//'   n = n,
//'   mu0 = mu0,
//'   sigma0_sqrt = sigma0_sqrt,
//'   alpha = alpha,
//'   beta = beta,
//'   psi_sqrt = psi_sqrt,
//'   nu = nu,
//'   lambda = lambda,
//'   theta_sqrt = theta_sqrt,
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
Rcpp::List SimSSM0Fixed(const int n, const arma::vec& mu0,
                        const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_sqrt,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_sqrt, const int time,
                        const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::mat y(num_manifest_vars, total_time);

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + sigma0_sqrt * arma::randn(num_latent_vars);
    y.col(0) =
        nu + lambda * eta.col(0) + theta_sqrt * arma::randn(num_manifest_vars);

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + beta * eta.col(t - 1) +
                   psi_sqrt * arma::randn(num_latent_vars);
      y.col(t) = nu + lambda * eta.col(t) +
                 theta_sqrt * arma::randn(num_manifest_vars);
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("time") = arma::regspace(1, time), Rcpp::Named("id") = id);
  }

  // Step 4: Return the results
  return out;
}
