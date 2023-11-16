// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-ou.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Data from the Ornstein–Uhlenbeck Model
//' using a State Space Model Parameterization (n = 1)
//'
//' This function simulates data from the Ornstein–Uhlenbeck model
//' using a state space model parameterization.
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
//'     \mathrm{d} \boldsymbol{\eta}_{t}
//'     =
//'     \boldsymbol{\Phi}
//'     \left(
//'     \boldsymbol{\mu}
//'     -
//'     \boldsymbol{\eta}_{t}
//'     \right)
//'     \mathrm{d}t
//'     +
//'     \boldsymbol{\Sigma}^{\frac{1}{2}}
//'     \mathrm{d}
//'     \mathbf{W}_{t}
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
//' @param mu Numeric vector.
//'   The long-term mean or equilibrium level
//'   (\eqn{\boldsymbol{\mu}}).
//' @param phi Numeric matrix.
//'   The rate of mean reversion,
//'   determining how quickly the variable returns to its mean
//'   (\eqn{\boldsymbol{\Phi}}).
//' @param sigma_sqrt Numeric matrix.
//'   Cholesky decomposition of the matrix of volatility
//'   or randomness in the process
//'   (\eqn{\boldsymbol{\Sigma}}).
//' @param delta_t Numeric.
//'   Time interval (\eqn{\delta_t}).
//' @inheritParams SimSSM0
//'
//' @references
//'   Chow, S.-M., Losardo, D., Park, J., & Molenaar, P. C. M. (2023).
//'   Continuous-time dynamic models: Connections to structural equation models and other discrete-time models.
//'   In R. H. Hoyle (Ed.), Handbook of structural equation modeling (2nd ed.). The Guilford Press.
//'
//'   Uhlenbeck, G. E., & Ornstein, L. S. (1930).
//'   On the theory of the brownian motion.
//'   *Physical Review*, *36*(5), 823–841.
//'   \doi{10.1103/physrev.36.823}
//'
//' @return Returns a list with the following elements:
//'   - `y`: A `t` by `k` matrix of values for the manifest variables.
//'   - `eta`: A `t` by `p` matrix of values for the latent variables.
//'   - `time`: A vector of continuous time points of length `t`
//'      starting from 0 with `delta_t` increments.
//'   - `n`: Number of individuals.
//'
//' @examples
//' # prepare parameters
//' set.seed(42)
//' p <- k <- 2
//' I <- diag(p)
//' I_sqrt <- chol(I)
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
//' ssm <- SimSSMOU(
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
//' @keywords simStateSpace sim
//' @export
// [[Rcpp::export]]
Rcpp::List SimSSMOU(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                    const arma::vec& mu, const arma::mat& phi,
                    const arma::mat& sigma_sqrt, const arma::vec& nu,
                    const arma::mat& lambda, const arma::mat& theta_sqrt,
                    const double delta_t, const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat y(num_manifest_vars, total_time);

  // Step 3: Get state space parameters
  arma::mat I = arma::eye<arma::mat>(num_latent_vars, num_latent_vars);
  arma::mat J = arma::eye<arma::mat>(num_latent_vars * num_latent_vars,
                                     num_latent_vars * num_latent_vars);
  arma::mat neg_phi = -1 * phi;
  // 3.1 beta
  arma::mat beta = arma::expmat(neg_phi * delta_t);  // A(Delta t)
  // 3.2 alpha
  arma::vec alpha = arma::inv(neg_phi) * (beta - I) * (phi * mu);  // b(Delta t)
  // 3.3 psi
  arma::mat neg_phi_hashtag = arma::kron(neg_phi, I) + arma::kron(I, neg_phi);
  arma::vec sigma_vec = arma::vectorise(sigma_sqrt * sigma_sqrt.t());
  arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                      (arma::expmat(neg_phi_hashtag * delta_t) - J) * sigma_vec;
  arma::mat psi_sqrt =
      arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

  // Step 4: Generate initial condition
  eta.col(0) = mu0 + sigma0_sqrt * arma::randn(num_latent_vars);
  y.col(0) =
      nu + lambda * eta.col(0) + theta_sqrt * arma::randn(num_manifest_vars);

  // Step 5: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) =
        alpha + beta * eta.col(t - 1) + psi_sqrt * arma::randn(num_latent_vars);
    y.col(t) =
        nu + lambda * eta.col(t) + theta_sqrt * arma::randn(num_manifest_vars);
  }

  // Step 6: If there is a burn-in period, remove it
  if (burn_in > 0) {
    y = y.cols(burn_in, total_time - 1);
    eta = eta.cols(burn_in, total_time - 1);
  }

  // Step 7: Return the transposed data matrices in a list
  return Rcpp::List::create(
      Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
      Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time));
}
