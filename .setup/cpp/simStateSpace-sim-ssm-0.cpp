// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Data from a State Space Model (n = 1)
//'
//' This function simulates data from a state space model.
//' See details for more information.
//'
//' @details The measurement model is given by
//'   \deqn{
//'     \mathbf{y}_{t}
//'     =
//'     \boldsymbol{\nu}
//'     +
//'     \boldsymbol{\Lambda}
//'     \boldsymbol{\eta}_{t}
//'     +
//'     \boldsymbol{\varepsilon}_{t}
//'     \quad
//'     \mathrm{with}
//'     \quad
//'     \boldsymbol{\varepsilon}_{t}
//'     \sim
//'     \mathcal{N}
//'     \left(
//'     \mathbf{0},
//'     \boldsymbol{\Theta}
//'     \right)
//'   }
//'   where \eqn{\mathbf{y}_{t}}, \eqn{\boldsymbol{\eta}_{t}},
//'   and \eqn{\boldsymbol{\varepsilon}_{t}}
//'   are random variables and \eqn{\boldsymbol{\nu}},
//'   \eqn{\boldsymbol{\Lambda}},
//'   and \eqn{\boldsymbol{\Theta}} are model parameters.
//'   \eqn{\mathbf{y}_{t}} is a vector of observed random variables
//'   at time \eqn{t},
//'   \eqn{\boldsymbol{\eta}_{t}} is a vector of latent random variables
//'   at time \eqn{t},
//'   and \eqn{\boldsymbol{\varepsilon}_{t}}
//'   is a vector of random measurement errors
//'   at time \eqn{t},
//'   while \eqn{\boldsymbol{\nu}} is a vector of intercept,
//'   \eqn{\boldsymbol{\Lambda}} is a matrix of factor loadings,
//'   and \eqn{\boldsymbol{\Theta}} is the covariance matrix of
//'   \eqn{\boldsymbol{\varepsilon}}.
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
//'   time \eqn{t - 1},
//'   and \eqn{\boldsymbol{\zeta}_{t}} is a vector of dynamic noise
//'   at time \eqn{t} while \eqn{\boldsymbol{\alpha}}
//'   is a vector of intercepts,
//'   \eqn{\boldsymbol{\beta}} is a matrix of autoregression
//'   and cross regression coefficients,
//'   and \eqn{\boldsymbol{\Psi}} is the covariance matrix of
//'   \eqn{\boldsymbol{\zeta}_{t}}.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param mu0 Numeric vector.
//'   Mean of initial latent variable values
//'   (\eqn{\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}}).
//' @param sigma0_sqrt Numeric matrix.
//'   Cholesky decomposition of the covariance matrix
//'   of initial latent variable values
//'   (\eqn{\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}}).
//' @param alpha Numeric vector.
//'   Vector of intercepts for the dynamic model
//'   (\eqn{\boldsymbol{\alpha}}).
//' @param beta Numeric matrix.
//'   Transition matrix relating the values of the latent variables
//'   at time `t - 1` to those at time `t`
//'   (\eqn{\boldsymbol{\beta}}).
//' @param psi_sqrt Numeric matrix.
//'   Cholesky decomposition of the process noise covariance matrix
//'   (\eqn{\boldsymbol{\Psi}}).
//' @param nu Numeric vector.
//'   Vector of intercepts for the measurement model
//'   (\eqn{\boldsymbol{\nu}}).
//' @param lambda Numeric matrix.
//'   Factor loading matrix linking the latent variables
//'   to the observed variables
//'   (\eqn{\boldsymbol{\Lambda}}).
//' @param theta_sqrt Numeric matrix.
//'   Cholesky decomposition of the measurement error covariance matrix
//'   (\eqn{\boldsymbol{\Theta}}).
//' @param time Positive integer.
//'   Number of time points to simulate.
//' @param burn_in Positive integer.
//'   Number of burn-in points to exclude before returning the results.
//'
//' @references
//'   Shumway, R. H., & Stoffer, D. S. (2017).
//'   *Time series analysis and its applications: With R examples*.
//'   Springer International Publishing.
//'   \doi{10.1007/978-3-319-52452-8}
//'
//' @return Returns a list with the following elements:
//'   - `y`: A `t` by `k` matrix of values for the manifest variables.
//'   - `eta`: A `t` by `p` matrix of values for the latent variables.
//'   - `time`: A vector of discrete time points from 1 to `t`.
//'   - `n`: Number of individuals.
//'
//' @examples
//' # prepare parameters
//' set.seed(42)
//' k <- p <- 3
//' I <- diag(k)
//' I_sqrt <- chol(I)
//' null_vec <- rep(x = 0, times = k)
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
//' ssm <- SimSSM0(
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
Rcpp::List SimSSM0(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                   const arma::vec& alpha, const arma::mat& beta,
                   const arma::mat& psi_sqrt, const arma::vec& nu,
                   const arma::mat& lambda, const arma::mat& theta_sqrt,
                   const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat y(num_manifest_vars, total_time);

  // Step 3: Generate initial condition
  eta.col(0) = mu0 + sigma0_sqrt * arma::randn(num_latent_vars);
  y.col(0) =
      nu + lambda * eta.col(0) + theta_sqrt * arma::randn(num_manifest_vars);

  // Step 4: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) =
        alpha + beta * eta.col(t - 1) + psi_sqrt * arma::randn(num_latent_vars);
    y.col(t) =
        nu + lambda * eta.col(t) + theta_sqrt * arma::randn(num_manifest_vars);
  }

  // Step 5: If there is a burn-in period, remove it
  if (burn_in > 0) {
    y = y.cols(burn_in, total_time - 1);
    eta = eta.cols(burn_in, total_time - 1);
  }

  // Step 6: Return the transposed data matrices in a list
  return Rcpp::List::create(Rcpp::Named("y") = y.t(),
                            Rcpp::Named("eta") = eta.t(),
                            Rcpp::Named("time") = arma::regspace(1, time));
}
