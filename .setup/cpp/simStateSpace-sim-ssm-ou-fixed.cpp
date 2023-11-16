// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-ou-fixed.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Data from an Ornstein–Uhlenbeck Model
//' using a State Space Model Parameterization
//' for n > 1 Individuals (Fixed Parameters)
//'
//' This function simulates data from an Ornstein–Uhlenbeck model
//' using a state space model parameterization
//' for `n > 1` individuals.
//' In this model,
//' the parameters are invariant across individuals.
//' See details for more information.
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
//'     \mathrm{d} \boldsymbol{\eta}_{i, t}
//'     =
//'     \boldsymbol{\Phi}
//'     \left(
//'     \boldsymbol{\mu}
//'     -
//'     \boldsymbol{\eta}_{i, t}
//'     \right)
//'     \mathrm{d}t
//'     +
//'     \boldsymbol{\Sigma}^{\frac{1}{2}}
//'     \mathrm{d}
//'     \mathbf{W}_{i, t}
//'   }
//'   where \eqn{\boldsymbol{\mu}} is the long-term mean or equilibrium level,
//'   \eqn{\boldsymbol{\Phi}} is the rate of mean reversion,
//'   determining how quickly the variable returns to its mean,
//'   \eqn{\boldsymbol{\Sigma}} is the matrix of volatility
//'   or randomness in the process, and \eqn{\mathrm{d}\boldsymbol{W}}
//'   is a Wiener process or Brownian motion,
//'   which represents random fluctuations.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @inheritParams SimSSMOU
//' @inheritParams SimSSM0Fixed
//' @inherit SimSSMOU references
//'
//' @return Returns a list of length `n`.
//'   Each element is a list with the following elements:
//'   - `y`: A `t` by `k` matrix of values for the manifest variables.
//'   - `eta`: A `t` by `p` matrix of values for the latent variables.
//'   - `time`: A vector of continuous time points of length `t`
//'      starting from 0 with `delta_t` increments.
//'   - `id`: A vector of ID numbers of length `t`.
//'   - `n`: Number of individuals.
//'
//' @examples
//' # prepare parameters
//' set.seed(42)
//' p <- k <- 2
//' I <- diag(p)
//' I_sqrt <- chol(I)
//' n <- 5
//' mu0 <- c(-3.0, 1.5)
//' sigma0_sqrt <- I_sqrt
//' mu <- c(5.76, 5.18)
//' phi <- matrix(data = c(0.10, -0.05, -0.05, 0.10), nrow = p)
//' sigma_sqrt <- chol(
//'   matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
//' )
//' nu <- rep(x = 0, times = k)
//' lambda <- diag(k)
//' theta_sqrt <- chol(diag(x = 0.50, nrow = k))
//' delta_t <- 0.10
//' time <- 50
//' burn_in <- 0
//'
//' # generate data
//' ssm <- SimSSMOUFixed(
//'   n = n,
//'   mu0 = mu0,
//'   sigma0_sqrt = sigma0_sqrt,
//'   mu = mu,
//'   phi = phi,
//'   sigma_sqrt = sigma_sqrt,
//'   nu = nu,
//'   lambda = lambda,
//'   theta_sqrt = theta_sqrt,
//'   delta_t = delta_t,
//'   time = time,
//'   burn_in = burn_in
//' )
//'
//' str(ssm)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim ou
//' @export
// [[Rcpp::export]]
Rcpp::List SimSSMOUFixed(const int n, const arma::vec& mu0,
                         const arma::mat& sigma0_sqrt, const arma::vec& mu,
                         const arma::mat& phi, const arma::mat& sigma_sqrt,
                         const arma::vec& nu, const arma::mat& lambda,
                         const arma::mat& theta_sqrt, const double delta_t,
                         const int time, const int burn_in) {
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

    // Step 3.2: Get state space parameters
    arma::mat I = arma::eye<arma::mat>(num_latent_vars, num_latent_vars);
    arma::mat J = arma::eye<arma::mat>(num_latent_vars * num_latent_vars,
                                       num_latent_vars * num_latent_vars);
    arma::mat neg_phi = -1 * phi;
    // 3.2.1 beta
    arma::mat beta = arma::expmat(neg_phi * delta_t);  // A(Delta t)
    // 3.2.2 alpha
    arma::vec alpha =
        arma::inv(neg_phi) * (beta - I) * (phi * mu);  // b(Delta t)

    // 3.2.3 psi
    arma::mat neg_phi_hashtag = arma::kron(neg_phi, I) + arma::kron(I, neg_phi);
    arma::vec sigma_vec = arma::vectorise(sigma_sqrt * sigma_sqrt.t());
    arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                        (arma::expmat(neg_phi_hashtag * delta_t) - J) *
                        sigma_vec;
    arma::mat psi_sqrt =
        arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0 + sigma0_sqrt * arma::randn(num_latent_vars);
    y.col(0) =
        nu + lambda * eta.col(0) + theta_sqrt * arma::randn(num_manifest_vars);

    // Step 3.4: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + beta * eta.col(t - 1) +
                   psi_sqrt * arma::randn(num_latent_vars);
      y.col(t) = nu + lambda * eta.col(t) +
                 theta_sqrt * arma::randn(num_manifest_vars);
    }

    // Step 3.5: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
    }

    // Step 3.6: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.7: Return the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
        Rcpp::Named("id") = id);
  }

  // Step 4: Return results
  return out;
}
