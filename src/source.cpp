// -----------------------------------------------------------------------------
// edit .setup/cpp/000-forward-declarations.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List SimSSM0(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                   const arma::vec& alpha, const arma::mat& beta,
                   const arma::mat& psi_sqrt, const arma::vec& nu,
                   const arma::mat& lambda, const arma::mat& theta_sqrt,
                   const int time, const int burn_in);

Rcpp::List SimSSMVAR(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                     const arma::vec& alpha, const arma::mat& beta,
                     const arma::mat& psi_sqrt, const int time,
                     const int burn_in);

Rcpp::List SimSSMOU(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                    const arma::vec& mu, const arma::mat& phi,
                    const arma::mat& sigma_sqrt, const arma::vec& nu,
                    const arma::mat& lambda, const arma::mat& theta_sqrt,
                    const double delta_t, const int time, const int burn_in);

Rcpp::List OU2SSM(const arma::vec& mu, const arma::mat& phi,
                  const arma::mat& sigma_sqrt, const double delta_t);

Rcpp::List SimSSM0Fixed(const int n, const arma::vec& mu0,
                        const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_sqrt,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_sqrt, const int time,
                        const int burn_in);

Rcpp::List SimSSMVARFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                          const arma::mat& beta, const arma::mat& psi_sqrt,
                          const int time, const int burn_in);

Rcpp::List SimSSMOUFixed(const int n, const arma::vec& mu0,
                         const arma::mat& sigma0_sqrt, const arma::vec& mu,
                         const arma::mat& phi, const arma::mat& sigma_sqrt,
                         const arma::vec& nu, const arma::mat& lambda,
                         const arma::mat& theta_sqrt, const double delta_t,
                         const int time, const int burn_in);
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-ou-2-ssm.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Convert Parameters from the Ornstein–Uhlenbeck Model
//' to State Space Model Parameterization
//'
//' This function converts parameters from the Ornstein–Uhlenbeck model
//' to state space model parameterization.
//' See details for more information.
//'
//' @details The state space parameters
//'   as a function of the  Ornstein–Uhlenbeck model parameters
//'   are given by
//'   \deqn{
//'       \boldsymbol{\beta}
//'       =
//'       \exp{
//'         \left(
//'           - \boldsymbol{\Phi}
//'           \Delta_{t}
//'         \right)
//'       }
//'   }
//'
//'   \deqn{
//'       \boldsymbol{\alpha}
//'       =
//'       - \boldsymbol{\Phi}^{-1}
//'       \left(
//'         \boldsymbol{\beta} - \mathbf{I}_{p}
//'       \right)
//'   }
//'
//'   \deqn{
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Psi}
//'       \right)
//'       =
//'       \left\{
//'         \left[
//'           \left(
//'             - \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
//'           \right)
//'           +
//'           \left(
//'             \mathbf{I}_{p} \otimes - \boldsymbol{\Phi}
//'           \right)
//'         \right]
//'         \left[
//'           \exp
//'           \left(
//'             \left[
//'               \left(
//'                 - \boldsymbol{\Phi} \otimes \mathbf{I}_{p}
//'               \right)
//'               +
//'               \left(
//'                 \mathbf{I}_{p} \otimes - \boldsymbol{\Phi}
//'               \right)
//'             \right]
//'             \Delta_{t}
//'         \right)
//'         -
//'         \mathbf{I}_{p \times p}
//'       \right]
//'       \mathrm{vec}
//'       \left(
//'         \boldsymbol{\Sigma}
//'       \right)
//'     \right\}
//'   }
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @inheritParams SimSSMOU
//'
//' @return Returns a list of state space parameters:
//'   - alpha: Numeric vector.
//'     Vector of intercepts for the dynamic model
//'     (\eqn{\boldsymbol{\alpha}}).
//'   - beta: Numeric matrix.
//'     Transition matrix relating the values of the latent variables
//'     at time `t - 1` to those at time `t`
//'     (\eqn{\boldsymbol{\beta}}).
//'   - psi: Numeric matrix.
//'     The process noise covariance matrix
//'     (\eqn{\boldsymbol{\Psi}}).
//'
//' @examples
//' p <- k <- 2
//' mu <- c(5.76, 5.18)
//' phi <- matrix(data = c(0.10, -0.05, -0.05, 0.10), nrow = p)
//' sigma_sqrt <- chol(
//'   matrix(data = c(2.79, 0.06, 0.06, 3.27), nrow = p)
//' )
//' delta_t <- 0.10
//'
//' OU2SSM(
//'   mu = mu,
//'   phi = phi,
//'   sigma_sqrt = sigma_sqrt,
//'   delta_t = delta_t
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim
//' @export
// [[Rcpp::export]]
Rcpp::List OU2SSM(const arma::vec& mu, const arma::mat& phi,
                  const arma::mat& sigma_sqrt, const double delta_t) {
  // Step 1: Determine indices
  int num_latent_vars = mu.n_elem;

  // Step 2: Get state space parameters
  arma::mat I = arma::eye<arma::mat>(num_latent_vars, num_latent_vars);
  arma::mat J = arma::eye<arma::mat>(num_latent_vars * num_latent_vars,
                                     num_latent_vars * num_latent_vars);
  arma::mat neg_phi = -1 * phi;
  // 2.1 beta
  arma::mat beta = arma::expmat(neg_phi * delta_t);  // A(Delta t)
  // 2.2 alpha
  arma::vec alpha = arma::inv(neg_phi) * (beta - I) * phi * mu;  // b(Delta t)
  // 2.3 psi
  arma::mat neg_phi_hashtag = arma::kron(neg_phi, I) + arma::kron(I, neg_phi);
  arma::vec sigma_vec = arma::vectorise(sigma_sqrt * sigma_sqrt.t());
  arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                      (arma::expmat(neg_phi_hashtag * delta_t) - J) * sigma_vec;
  arma::mat psi = arma::reshape(psi_vec, num_latent_vars, num_latent_vars);

  // Step 3: Return state space parameters in a list
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi") = psi);
}
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
//' @keywords simStateSpace sim
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
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-var-fixed.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Data from a Vector Autoregressive Model
//' using a State Space Model Parameterization
//' for n > 1 Individuals (Fixed Parameters)
//'
//' This function simulates data from a vector autoregressive model
//' using a state space model parameterization
//' for `n > 1` individuals.
//' In this model,
//' the parameters are invariant across individuals.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @inheritParams SimSSM0Fixed
//' @inherit SimSSM0Fixed return
//' @inherit SimSSM0 references
//'
//' @examples
//' # prepare parameters
//' set.seed(42)
//' k <- 3
//' iden <- diag(k)
//' iden_sqrt <- chol(iden)
//' null_vec <- rep(x = 0, times = k)
//' n <- 5
//' mu0 <- null_vec
//' sigma0_sqrt <- iden_sqrt
//' alpha <- null_vec
//' beta <- diag(x = 0.5, nrow = k)
//' psi_sqrt <- iden_sqrt
//' time <- 50
//' burn_in <- 0
//'
//' ssm <- SimSSMVARFixed(
//'   n = n,
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
Rcpp::List SimSSMVARFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                          const arma::mat& beta, const arma::mat& psi_sqrt,
                          const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + sigma0_sqrt * arma::randn(num_latent_vars);

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + beta * eta.col(t - 1) +
                   psi_sqrt * arma::randn(num_latent_vars);
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      eta = eta.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Return the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("y") = eta.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("time") = arma::regspace(1, time), Rcpp::Named("id") = id);
  }

  // Step 4: Return the results
  return out;
}
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
