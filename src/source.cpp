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
//'   - `alpha`: Numeric vector.
//'     Vector of intercepts for the dynamic model
//'     (\eqn{\boldsymbol{\alpha}}).
//'   - `beta`: Numeric matrix.
//'     Transition matrix relating the values of the latent variables
//'     at time `t - 1` to those at time `t`
//'     (\eqn{\boldsymbol{\beta}}).
//'   - `psi`: Numeric matrix.
//'     The process noise covariance matrix
//'     (\eqn{\boldsymbol{\Psi}}).
//'
//' @examples
//' p <- k <- 2
//' mu <- c(5.76, 5.18)
//' phi <- matrix(
//'   data = c(0.10, -0.05, -0.05, 0.10),
//'   nrow = p
//' )
//' sigma <- matrix(
//'   data = c(2.79, 0.06, 0.06, 3.27),
//'   nrow = p
//' )
//' delta_t <- 0.10
//'
//' OU2SSM(
//'   mu = mu,
//'   phi = phi,
//'   sigma = sigma,
//'   delta_t = delta_t
//' )
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim ou
//' @export
// [[Rcpp::export]]
Rcpp::List OU2SSM(const arma::vec& mu, const arma::mat& phi,
                  const arma::mat& sigma, const double delta_t) {
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
  arma::vec sigma_vec = arma::vectorise(sigma);
  arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                      (arma::expmat(neg_phi_hashtag * delta_t) - J) * sigma_vec;
  arma::mat psi = arma::reshape(psi_vec, num_latent_vars, num_latent_vars);

  // Step 3: Return state space parameters in a list
  return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                            Rcpp::Named("beta") = beta,
                            Rcpp::Named("psi") = psi);
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0)]]
Rcpp::List SimSSM0(const arma::vec& mu0, const arma::mat& sigma0_l,
                   const arma::vec& alpha, const arma::mat& beta,
                   const arma::mat& psi_l, const arma::vec& nu,
                   const arma::mat& lambda, const arma::mat& theta_l,
                   const int time, const int burn_in) {
  // Note:
  // sigma0_l, psi_l, and theta_l are L in A = L * L^T
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat y(num_manifest_vars, total_time);
  arma::vec id(total_time, arma::fill::ones);

  // Step 3: Generate initial condition
  eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars));
  y.col(0) =
      nu + (lambda * eta.col(0)) + (theta_l * arma::randn(num_manifest_vars));

  // Step 4: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                 (psi_l * arma::randn(num_latent_vars));
    y.col(t) =
        nu + (lambda * eta.col(t)) + (theta_l * arma::randn(num_manifest_vars));
  }

  // Step 5: If there is a burn-in period, remove it
  if (burn_in > 0) {
    y = y.cols(burn_in, total_time - 1);
    eta = eta.cols(burn_in, total_time - 1);
    id = id.subvec(burn_in, total_time - 1);
  }

  // Step 6: Return the transposed data matrices in a list
  return Rcpp::List::create(
      Rcpp::Named("id") = id, Rcpp::Named("time") = arma::regspace(0, time - 1),
      Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-fixed-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0Fixed)]]
Rcpp::List SimSSM0Fixed(const int n, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l, const int time,
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
    eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars));
    y.col(0) =
        nu + (lambda * eta.col(0)) + (theta_l * arma::randn(num_manifest_vars));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars));
      y.col(t) = nu + (lambda * eta.col(t)) +
                 (theta_l * arma::randn(num_manifest_vars));
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
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0IVary)]]
Rcpp::List SimSSM0IVary(const int n, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l, const int time,
                        const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  arma::vec mu0_temp = mu0[0];
  arma::vec nu_temp = nu[0];
  int num_latent_vars = mu0_temp.n_elem;
  int num_manifest_vars = nu_temp.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::mat y(num_manifest_vars, total_time);
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    arma::vec alpha_temp = alpha[i];
    arma::mat beta_temp = beta[i];
    arma::mat psi_l_temp = psi_l[i];
    arma::vec nu_temp = nu[i];
    arma::mat lambda_temp = lambda[i];
    arma::mat theta_l_temp = theta_l[i];

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(num_latent_vars));
    y.col(0) = nu_temp + (lambda_temp * eta.col(0)) +
               (theta_l_temp * arma::randn(num_manifest_vars));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha_temp + (beta_temp * eta.col(t - 1)) +
                   (psi_l_temp * arma::randn(num_latent_vars));
      y.col(t) = nu_temp + (lambda_temp * eta.col(t)) +
                 (theta_l_temp * arma::randn(num_manifest_vars));
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
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-lin-growth-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0LinGrowth)]]
Rcpp::List SimSSM0LinGrowth(const int n, const arma::vec& mu0,
                            const arma::mat& sigma0_l, const double theta_l,
                            const int time) {
  // Step 1: Create constant vectors and matrices
  arma::mat lambda = {{1, 0}};
  arma::mat beta = {{1, 1}, {0, 1}};

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(2, time);
    arma::mat y(1, time);

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(2));
    y.col(0) = (lambda * eta.col(0)) + (theta_l * arma::randn(1));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = (beta * eta.col(t - 1));
      y.col(t) = (lambda * eta.col(t)) + (theta_l * arma::randn(1));
    }

    // Step 3.4: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.5: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-lin-growth-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0LinGrowthIVary)]]
Rcpp::List SimSSM0LinGrowthIVary(const int n, const Rcpp::List& mu0,
                                 const Rcpp::List& sigma0_l,
                                 const Rcpp::List& theta_l, const int time) {
  // Step 1: Create constant vectors and matrices
  arma::mat lambda = {{1, 0}};
  arma::mat beta = {{1, 1}, {0, 1}};

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(2, time);
    arma::mat y(1, time);
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    double theta_l_temp = theta_l[i];

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(2));
    y.col(0) = (lambda * eta.col(0)) + (theta_l_temp * arma::randn(1));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = (beta * eta.col(t - 1));
      y.col(t) = (lambda * eta.col(t)) + (theta_l_temp * arma::randn(1));
    }

    // Step 3.4: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.5: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-ou-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0OU)]]
Rcpp::List SimSSM0OU(const arma::vec& mu0, const arma::mat& sigma0_l,
                     const arma::vec& mu, const arma::mat& phi,
                     const arma::mat& sigma_l, const arma::vec& nu,
                     const arma::mat& lambda, const arma::mat& theta_l,
                     const double delta_t, const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat y(num_manifest_vars, total_time);
  arma::vec id(total_time, arma::fill::ones);

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
  arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
  arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                      (arma::expmat(neg_phi_hashtag * delta_t) - J) * sigma_vec;
  arma::mat psi_l =
      arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

  // Step 4: Generate initial condition
  eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars));
  y.col(0) =
      nu + (lambda * eta.col(0)) + (theta_l * arma::randn(num_manifest_vars));

  // Step 5: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                 (psi_l * arma::randn(num_latent_vars));
    y.col(t) =
        nu + (lambda * eta.col(t)) + (theta_l * arma::randn(num_manifest_vars));
  }

  // Step 6: If there is a burn-in period, remove it
  if (burn_in > 0) {
    y = y.cols(burn_in, total_time - 1);
    eta = eta.cols(burn_in, total_time - 1);
    id = id.subvec(burn_in, total_time - 1);
  }

  // Step 7: Return the transposed data matrices in a list
  return Rcpp::List::create(
      Rcpp::Named("id") = id,
      Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
      Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-ou-fixed-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0OUFixed)]]
Rcpp::List SimSSM0OUFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_l, const arma::vec& mu,
                          const arma::mat& phi, const arma::mat& sigma_l,
                          const arma::vec& nu, const arma::mat& lambda,
                          const arma::mat& theta_l, const double delta_t,
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
    arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
    arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                        (arma::expmat(neg_phi_hashtag * delta_t) - J) *
                        sigma_vec;
    arma::mat psi_l =
        arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars));
    y.col(0) =
        nu + (lambda * eta.col(0)) + (theta_l * arma::randn(num_manifest_vars));

    // Step 3.4: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars));
      y.col(t) = nu + (lambda * eta.col(t)) +
                 (theta_l * arma::randn(num_manifest_vars));
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
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-ou-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0OUIVary)]]
Rcpp::List SimSSM0OUIVary(const int n, const Rcpp::List& mu0,
                          const Rcpp::List& sigma0_l, const Rcpp::List& mu,
                          const Rcpp::List& phi, const Rcpp::List& sigma_l,
                          const Rcpp::List& nu, const Rcpp::List& lambda,
                          const Rcpp::List& theta_l, const double delta_t,
                          const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  arma::vec mu0_temp = mu0[0];
  arma::vec nu_temp = nu[0];
  int num_latent_vars = mu0_temp.n_elem;
  int num_manifest_vars = nu_temp.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::mat y(num_manifest_vars, total_time);
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    arma::vec mu_temp = mu[i];
    arma::mat phi_temp = phi[i];
    arma::mat sigma_l_temp = sigma_l[i];
    arma::vec nu_temp = nu[i];
    arma::mat lambda_temp = lambda[i];
    arma::mat theta_l_temp = theta_l[i];

    // Step 3.2: Get state space parameters
    arma::mat I = arma::eye<arma::mat>(num_latent_vars, num_latent_vars);
    arma::mat J = arma::eye<arma::mat>(num_latent_vars * num_latent_vars,
                                       num_latent_vars * num_latent_vars);
    arma::mat neg_phi = -1 * phi_temp;
    // 3.2.1 beta_temp
    arma::mat beta_temp = arma::expmat(neg_phi * delta_t);  // A(Delta t)
    // 3.2.2 alpha_temp
    arma::vec alpha_temp = arma::inv(neg_phi) * (beta_temp - I) *
                           (phi_temp * mu_temp);  // b(Delta t)

    // 3.2.3 psi
    arma::mat neg_phi_hashtag = arma::kron(neg_phi, I) + arma::kron(I, neg_phi);
    arma::vec sigma_vec = arma::vectorise(sigma_l_temp * sigma_l_temp.t());
    arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                        (arma::expmat(neg_phi_hashtag * delta_t) - J) *
                        sigma_vec;
    arma::mat psi_l =
        arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(num_latent_vars));
    y.col(0) = nu_temp + (lambda_temp * eta.col(0)) +
               (theta_l_temp * arma::randn(num_manifest_vars));

    // Step 3.4: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha_temp + (beta_temp * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars));
      y.col(t) = nu_temp + (lambda_temp * eta.col(t)) +
                 (theta_l_temp * arma::randn(num_manifest_vars));
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
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-var-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0VAR)]]
Rcpp::List SimSSM0VAR(const arma::vec& mu0, const arma::mat& sigma0_l,
                      const arma::vec& alpha, const arma::mat& beta,
                      const arma::mat& psi_l, const int time,
                      const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::vec id(total_time, arma::fill::ones);

  // Step 3: Generate initial condition
  eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars));

  // Step 4: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                 (psi_l * arma::randn(num_latent_vars));
  }

  // Step 5: If there is a burn-in period, remove it
  if (burn_in > 0) {
    eta = eta.cols(burn_in, total_time - 1);
    id = id.subvec(burn_in, total_time - 1);
  }

  // Step 6: Return the transposed data matrices in a list
  return Rcpp::List::create(
      Rcpp::Named("id") = id, Rcpp::Named("time") = arma::regspace(0, time - 1),
      Rcpp::Named("y") = eta.t(), Rcpp::Named("eta") = eta.t());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-var-fixed-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0VARFixed)]]
Rcpp::List SimSSM0VARFixed(const int n, const arma::vec& mu0,
                           const arma::mat& sigma0_l, const arma::vec& alpha,
                           const arma::mat& beta, const arma::mat& psi_l,
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
    eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars));
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
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = eta.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-var-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0VARIVary)]]
Rcpp::List SimSSM0VARIVary(const int n, const Rcpp::List& mu0,
                           const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                           const Rcpp::List& beta, const Rcpp::List& psi_l,
                           const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  arma::vec mu0_temp = mu0[0];
  int num_latent_vars = mu0_temp.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    arma::vec alpha_temp = alpha[i];
    arma::mat beta_temp = beta[i];
    arma::mat psi_l_temp = psi_l[i];

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(num_latent_vars));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha_temp + (beta_temp * eta.col(t - 1)) +
                   psi_l_temp * arma::randn(num_latent_vars);
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      eta = eta.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = eta.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1)]]
Rcpp::List SimSSM1(const arma::vec& mu0, const arma::mat& sigma0_l,
                   const arma::vec& alpha, const arma::mat& beta,
                   const arma::mat& psi_l, const arma::vec& nu,
                   const arma::mat& lambda, const arma::mat& theta_l,
                   const arma::mat& gamma_eta, const arma::mat& x,
                   const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat y(num_manifest_vars, total_time);
  arma::mat x_t = x.t();
  arma::vec id(total_time, arma::fill::ones);

  // Step 3: Generate initial condition
  eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
               (gamma_eta * x_t.col(0));
  y.col(0) =
      nu + (lambda * eta.col(0)) + (theta_l * arma::randn(num_manifest_vars));

  // Step 4: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                 (psi_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(t));
    y.col(t) =
        nu + (lambda * eta.col(t)) + (theta_l * arma::randn(num_manifest_vars));
  }

  // Step 5: If there is a burn-in period, remove it
  if (burn_in > 0) {
    y = y.cols(burn_in, total_time - 1);
    eta = eta.cols(burn_in, total_time - 1);
    x_t = x_t.cols(burn_in, total_time - 1);
    id = id.subvec(burn_in, total_time - 1);
  }

  // Step 6: Return the transposed data matrices in a list
  return Rcpp::List::create(
      Rcpp::Named("id") = id, Rcpp::Named("time") = arma::regspace(0, time - 1),
      Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
      Rcpp::Named("x") = x_t.t());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-fixed-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1Fixed)]]
Rcpp::List SimSSM1Fixed(const int n, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l, const arma::mat& gamma_eta,
                        const Rcpp::List& x, const int time,
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
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(0));
    y.col(0) =
        nu + (lambda * eta.col(0)) + (theta_l * arma::randn(num_manifest_vars));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars)) +
                   (gamma_eta * x_t.col(t));
      y.col(t) = nu + (lambda * eta.col(t)) +
                 (theta_l * arma::randn(num_manifest_vars));
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1IVary)]]
Rcpp::List SimSSM1IVary(const int n, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l, const Rcpp::List& gamma_eta,
                        const Rcpp::List& x, const int time,
                        const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  arma::vec mu0_temp = mu0[0];
  arma::vec nu_temp = nu[0];
  int num_latent_vars = mu0_temp.n_elem;
  int num_manifest_vars = nu_temp.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::mat y(num_manifest_vars, total_time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    arma::vec alpha_temp = alpha[i];
    arma::mat beta_temp = beta[i];
    arma::mat psi_l_temp = psi_l[i];
    arma::vec nu_temp = nu[i];
    arma::mat lambda_temp = lambda[i];
    arma::mat theta_l_temp = theta_l[i];
    arma::mat gamma_eta_temp = gamma_eta[i];

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(num_latent_vars)) +
                 (gamma_eta_temp * x_t.col(0));
    y.col(0) = nu_temp + (lambda_temp * eta.col(0)) +
               (theta_l_temp * arma::randn(num_manifest_vars));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha_temp + (beta_temp * eta.col(t - 1)) +
                   psi_l_temp * arma::randn(num_latent_vars) +
                   (gamma_eta_temp * x_t.col(t));
      y.col(t) = nu_temp + (lambda_temp * eta.col(t)) +
                 (theta_l_temp * arma::randn(num_manifest_vars));
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-lin-growth-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1LinGrowth)]]
Rcpp::List SimSSM1LinGrowth(const int n, const arma::vec& mu0,
                            const arma::mat& sigma0_l, const double theta_l,
                            const arma::mat& gamma_eta, const Rcpp::List& x,
                            const int time) {
  // Step 1: Create constant vectors and matrices
  arma::mat lambda = {{1, 0}};
  arma::mat beta = {{1, 1}, {0, 1}};

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(2, time);
    arma::mat y(1, time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(2)) + (gamma_eta * x_t.col(0));
    y.col(0) = (lambda * eta.col(0)) + (theta_l * arma::randn(1));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = (beta * eta.col(t - 1)) + (gamma_eta * x_t.col(t));
      y.col(t) = (lambda * eta.col(t)) + (theta_l * arma::randn(1));
    }

    // Step 3.4: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.5: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-lin-growth-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1LinGrowthIVary)]]
Rcpp::List SimSSM1LinGrowthIVary(const int n, const Rcpp::List& mu0,
                                 const Rcpp::List& sigma0_l,
                                 const Rcpp::List& theta_l,
                                 const Rcpp::List& gamma_eta,
                                 const Rcpp::List& x, const int time) {
  // Step 1: Create constant vectors and matrices
  arma::mat lambda = {{1, 0}};
  arma::mat beta = {{1, 1}, {0, 1}};

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(2, time);
    arma::mat y(1, time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    double theta_l_temp = theta_l[i];
    arma::mat gamma_eta_temp = gamma_eta[i];

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(2)) +
                 (gamma_eta_temp * x_t.col(0));
    y.col(0) = (lambda * eta.col(0)) + (theta_l_temp * arma::randn(1));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = (beta * eta.col(t - 1)) + (gamma_eta_temp * x_t.col(t));
      y.col(t) = (lambda * eta.col(t)) + (theta_l_temp * arma::randn(1));
    }

    // Step 3.4: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.5: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-ou-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1OU)]]
Rcpp::List SimSSM1OU(const arma::vec& mu0, const arma::mat& sigma0_l,
                     const arma::vec& mu, const arma::mat& phi,
                     const arma::mat& sigma_l, const arma::vec& nu,
                     const arma::mat& lambda, const arma::mat& theta_l,
                     const arma::mat& gamma_eta, const arma::mat& x,
                     const double delta_t, const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat y(num_manifest_vars, total_time);
  arma::mat x_t = x.t();
  arma::vec id(total_time, arma::fill::ones);

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
  arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
  arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                      (arma::expmat(neg_phi_hashtag * delta_t) - J) * sigma_vec;
  arma::mat psi_l =
      arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

  // Step 4: Generate initial condition
  eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
               (gamma_eta * x_t.col(0));
  y.col(0) =
      nu + (lambda * eta.col(0)) + (theta_l * arma::randn(num_manifest_vars));

  // Step 5: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                 (psi_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(t));
    y.col(t) =
        nu + (lambda * eta.col(t)) + (theta_l * arma::randn(num_manifest_vars));
  }

  // Step 6: If there is a burn-in period, remove it
  if (burn_in > 0) {
    y = y.cols(burn_in, total_time - 1);
    eta = eta.cols(burn_in, total_time - 1);
    x_t = x_t.cols(burn_in, total_time - 1);
    id = id.subvec(burn_in, total_time - 1);
  }

  // Step 7: Return the transposed data matrices in a list
  return Rcpp::List::create(
      Rcpp::Named("id") = id,
      Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
      Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
      Rcpp::Named("x") = x_t.t());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-ou-fixed-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1OUFixed)]]
Rcpp::List SimSSM1OUFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_l, const arma::vec& mu,
                          const arma::mat& phi, const arma::mat& sigma_l,
                          const arma::vec& nu, const arma::mat& lambda,
                          const arma::mat& theta_l, const arma::mat& gamma_eta,
                          const Rcpp::List& x, const double delta_t,
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
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();

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
    arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
    arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                        (arma::expmat(neg_phi_hashtag * delta_t) - J) *
                        sigma_vec;
    arma::mat psi_l =
        arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(0));
    y.col(0) =
        nu + (lambda * eta.col(0)) + (theta_l * arma::randn(num_manifest_vars));

    // Step 3.4: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars)) +
                   (gamma_eta * x_t.col(t));
      y.col(t) = nu + (lambda * eta.col(t)) +
                 (theta_l * arma::randn(num_manifest_vars));
    }

    // Step 3.5: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.6: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.7: Return the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-ou-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1OUIVary)]]
Rcpp::List SimSSM1OUIVary(const int n, const Rcpp::List& mu0,
                          const Rcpp::List& sigma0_l, const Rcpp::List& mu,
                          const Rcpp::List& phi, const Rcpp::List& sigma_l,
                          const Rcpp::List& nu, const Rcpp::List& lambda,
                          const Rcpp::List& theta_l,
                          const Rcpp::List& gamma_eta, const Rcpp::List& x,
                          const double delta_t, const int time,
                          const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  arma::vec mu0_temp = mu0[0];
  arma::vec nu_temp = nu[0];
  int num_latent_vars = mu0_temp.n_elem;
  int num_manifest_vars = nu_temp.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::mat y(num_manifest_vars, total_time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    arma::vec mu_temp = mu[i];
    arma::mat phi_temp = phi[i];
    arma::mat sigma_l_temp = sigma_l[i];
    arma::vec nu_temp = nu[i];
    arma::mat lambda_temp = lambda[i];
    arma::mat theta_l_temp = theta_l[i];
    arma::mat gamma_eta_temp = gamma_eta[i];

    // Step 3.2: Get state space parameters
    arma::mat I = arma::eye<arma::mat>(num_latent_vars, num_latent_vars);
    arma::mat J = arma::eye<arma::mat>(num_latent_vars * num_latent_vars,
                                       num_latent_vars * num_latent_vars);
    arma::mat neg_phi = -1 * phi_temp;
    // 3.2.1 beta_temp
    arma::mat beta_temp = arma::expmat(neg_phi * delta_t);  // A(Delta t)
    // 3.2.2 alpha_temp
    arma::vec alpha_temp = arma::inv(neg_phi) * (beta_temp - I) *
                           (phi_temp * mu_temp);  // b(Delta t)

    // 3.2.3 psi
    arma::mat neg_phi_hashtag = arma::kron(neg_phi, I) + arma::kron(I, neg_phi);
    arma::vec sigma_vec = arma::vectorise(sigma_l_temp * sigma_l_temp.t());
    arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                        (arma::expmat(neg_phi_hashtag * delta_t) - J) *
                        sigma_vec;
    arma::mat psi_l =
        arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(num_latent_vars)) +
                 (gamma_eta_temp * x_t.col(0));
    y.col(0) = nu_temp + (lambda_temp * eta.col(0)) +
               (theta_l_temp * arma::randn(num_manifest_vars));

    // Step 3.4: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha_temp + (beta_temp * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars)) +
                   (gamma_eta_temp * x_t.col(t));
      y.col(t) = nu_temp + (lambda_temp * eta.col(t)) +
                 (theta_l_temp * arma::randn(num_manifest_vars));
    }

    // Step 3.5: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.6: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.7: Return the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-var-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1VAR)]]
Rcpp::List SimSSM1VAR(const arma::vec& mu0, const arma::mat& sigma0_l,
                      const arma::vec& alpha, const arma::mat& beta,
                      const arma::mat& psi_l, const arma::mat& gamma_eta,
                      const arma::mat& x, const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat x_t = x.t();
  arma::vec id(total_time, arma::fill::ones);

  // Step 3: Generate initial condition
  eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
               (gamma_eta * x_t.col(0));

  // Step 4: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                 (psi_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(t));
  }

  // Step 5: If there is a burn-in period, remove it
  if (burn_in > 0) {
    eta = eta.cols(burn_in, total_time - 1);
    x_t = x_t.cols(burn_in, total_time - 1);
    id = id.subvec(burn_in, total_time - 1);
  }

  // Step 6: Return the transposed data matrices in a list
  return Rcpp::List::create(
      Rcpp::Named("id") = id, Rcpp::Named("time") = arma::regspace(0, time - 1),
      Rcpp::Named("y") = eta.t(), Rcpp::Named("eta") = eta.t(),
      Rcpp::Named("x") = x_t.t());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-var-fixed-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1VARFixed)]]
Rcpp::List SimSSM1VARFixed(const int n, const arma::vec& mu0,
                           const arma::mat& sigma0_l, const arma::vec& alpha,
                           const arma::mat& beta, const arma::mat& psi_l,
                           arma::mat& gamma_eta, const Rcpp::List& x,
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
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(0));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars)) +
                   (gamma_eta * x_t.col(t));
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Return the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = eta.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-var-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1VARIVary)]]
Rcpp::List SimSSM1VARIVary(const int n, const Rcpp::List& mu0,
                           const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                           const Rcpp::List& beta, const Rcpp::List& psi_l,
                           const Rcpp::List& gamma_eta, const Rcpp::List& x,
                           const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  arma::vec mu0_temp = mu0[0];
  int num_latent_vars = mu0_temp.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    arma::vec alpha_temp = alpha[i];
    arma::mat beta_temp = beta[i];
    arma::mat psi_l_temp = psi_l[i];
    arma::mat gamma_eta_temp = gamma_eta[i];

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(num_latent_vars)) +
                 (gamma_eta_temp * x_t.col(0));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha_temp + (beta_temp * eta.col(t - 1)) +
                   psi_l_temp * arma::randn(num_latent_vars) +
                   (gamma_eta_temp * x_t.col(t));
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = eta.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-2-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM2)]]
Rcpp::List SimSSM2(const arma::vec& mu0, const arma::mat& sigma0_l,
                   const arma::vec& alpha, const arma::mat& beta,
                   const arma::mat& psi_l, const arma::vec& nu,
                   const arma::mat& lambda, const arma::mat& theta_l,
                   const arma::mat& gamma_y, const arma::mat& gamma_eta,
                   const arma::mat& x, const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat y(num_manifest_vars, total_time);
  arma::mat x_t = x.t();
  arma::vec id(total_time, arma::fill::ones);

  // Step 3: Generate initial condition
  eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
               (gamma_eta * x_t.col(0));
  y.col(0) = nu + (lambda * eta.col(0)) +
             (theta_l * arma::randn(num_manifest_vars)) +
             (gamma_y * x_t.col(0));

  // Step 4: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                 (psi_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(t));
    y.col(t) = nu + (lambda * eta.col(t)) +
               (theta_l * arma::randn(num_manifest_vars)) +
               (gamma_y * x_t.col(t));
  }

  // Step 5: If there is a burn-in period, remove it
  if (burn_in > 0) {
    y = y.cols(burn_in, total_time - 1);
    eta = eta.cols(burn_in, total_time - 1);
    x_t = x_t.cols(burn_in, total_time - 1);
    id = id.subvec(burn_in, total_time - 1);
  }

  // Step 6: Return the transposed data matrices in a list
  return Rcpp::List::create(
      Rcpp::Named("id") = id, Rcpp::Named("time") = arma::regspace(0, time - 1),
      Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
      Rcpp::Named("x") = x_t.t());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-2-fixed-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM2Fixed)]]
Rcpp::List SimSSM2Fixed(const int n, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l, const arma::mat& gamma_y,
                        const arma::mat& gamma_eta, const Rcpp::List& x,
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
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(0));
    y.col(0) = nu + (lambda * eta.col(0)) +
               (theta_l * arma::randn(num_manifest_vars)) +
               (gamma_y * x_t.col(0));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars)) +
                   (gamma_eta * x_t.col(t));
      y.col(t) = nu + (lambda * eta.col(t)) +
                 (theta_l * arma::randn(num_manifest_vars)) +
                 (gamma_y * x_t.col(t));
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-2-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM2IVary)]]
Rcpp::List SimSSM2IVary(const int n, const Rcpp::List& mu0,
                        const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                        const Rcpp::List& beta, const Rcpp::List& psi_l,
                        const Rcpp::List& nu, const Rcpp::List& lambda,
                        const Rcpp::List& theta_l, const Rcpp::List& gamma_y,
                        const Rcpp::List& gamma_eta, const Rcpp::List& x,
                        const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  arma::vec mu0_temp = mu0[0];
  arma::vec nu_temp = nu[0];
  int num_latent_vars = mu0_temp.n_elem;
  int num_manifest_vars = nu_temp.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::mat y(num_manifest_vars, total_time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    arma::vec alpha_temp = alpha[i];
    arma::mat beta_temp = beta[i];
    arma::mat psi_l_temp = psi_l[i];
    arma::vec nu_temp = nu[i];
    arma::mat lambda_temp = lambda[i];
    arma::mat theta_l_temp = theta_l[i];
    arma::mat gamma_y_temp = gamma_y[i];
    arma::mat gamma_eta_temp = gamma_eta[i];

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(num_latent_vars)) +
                 (gamma_eta_temp * x_t.col(0));
    y.col(0) = nu_temp + (lambda_temp * eta.col(0)) +
               (theta_l_temp * arma::randn(num_manifest_vars)) +
               (gamma_y_temp * x_t.col(0));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha_temp + (beta_temp * eta.col(t - 1)) +
                   psi_l_temp * arma::randn(num_latent_vars) +
                   (gamma_eta_temp * x_t.col(t));
      y.col(t) = nu_temp + (lambda_temp * eta.col(t)) +
                 (theta_l_temp * arma::randn(num_manifest_vars)) +
                 (gamma_y_temp * x_t.col(t));
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-2-lin-growth-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM2LinGrowth)]]
Rcpp::List SimSSM2LinGrowth(const int n, const arma::vec& mu0,
                            const arma::mat& sigma0_l, const double theta_l,
                            const arma::mat& gamma_y,
                            const arma::mat& gamma_eta, const Rcpp::List& x,
                            const int time) {
  // Step 1: Create constant vectors and matrices
  arma::mat lambda = {{1, 0}};
  arma::mat beta = {{1, 1}, {0, 1}};

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(2, time);
    arma::mat y(1, time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(2)) + (gamma_eta * x_t.col(0));
    y.col(0) = (lambda * eta.col(0)) + (theta_l * arma::randn(1)) +
               (gamma_y * x_t.col(0));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = (beta * eta.col(t - 1)) + (gamma_eta * x_t.col(t));
      y.col(t) = (lambda * eta.col(t)) + (theta_l * arma::randn(1)) +
                 (gamma_y * x_t.col(t));
    }

    // Step 3.4: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.5: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-2-lin-growth-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM2LinGrowthIVary)]]
Rcpp::List SimSSM2LinGrowthIVary(const int n, const Rcpp::List& mu0,
                                 const Rcpp::List& sigma0_l,
                                 const Rcpp::List& theta_l,
                                 const Rcpp::List& gamma_y,
                                 const Rcpp::List& gamma_eta,
                                 const Rcpp::List& x, const int time) {
  // Step 1: Create constant vectors and matrices
  arma::mat lambda = {{1, 0}};
  arma::mat beta = {{1, 1}, {0, 1}};

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(2, time);
    arma::mat y(1, time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    double theta_l_temp = theta_l[i];
    arma::mat gamma_y_temp = gamma_y[i];
    arma::mat gamma_eta_temp = gamma_eta[i];

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(2)) +
                 (gamma_eta_temp * x_t.col(0));
    y.col(0) = (lambda * eta.col(0)) + (theta_l_temp * arma::randn(1)) +
               (gamma_y_temp * x_t.col(0));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = (beta * eta.col(t - 1)) + (gamma_eta_temp * x_t.col(t));
      y.col(t) = (lambda * eta.col(t)) + (theta_l_temp * arma::randn(1)) +
                 (gamma_y_temp * x_t.col(t));
    }

    // Step 3.4: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.5: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return the results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-2-ou-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM2OU)]]
Rcpp::List SimSSM2OU(const arma::vec& mu0, const arma::mat& sigma0_l,
                     const arma::vec& mu, const arma::mat& phi,
                     const arma::mat& sigma_l, const arma::vec& nu,
                     const arma::mat& lambda, const arma::mat& theta_l,
                     const arma::mat& gamma_y, const arma::mat& gamma_eta,
                     const arma::mat& x, const double delta_t, const int time,
                     const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat y(num_manifest_vars, total_time);
  arma::mat x_t = x.t();
  arma::vec id(total_time, arma::fill::ones);

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
  arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
  arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                      (arma::expmat(neg_phi_hashtag * delta_t) - J) * sigma_vec;
  arma::mat psi_l =
      arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

  // Step 4: Generate initial condition
  eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
               (gamma_eta * x_t.col(0));
  y.col(0) = nu + (lambda * eta.col(0)) +
             (theta_l * arma::randn(num_manifest_vars)) +
             (gamma_y * x_t.col(0));

  // Step 5: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                 (psi_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(t));
    y.col(t) = nu + (lambda * eta.col(t)) +
               (theta_l * arma::randn(num_manifest_vars)) +
               (gamma_y * x_t.col(t));
  }

  // Step 6: If there is a burn-in period, remove it
  if (burn_in > 0) {
    y = y.cols(burn_in, total_time - 1);
    eta = eta.cols(burn_in, total_time - 1);
    x_t = x_t.cols(burn_in, total_time - 1);
    id = id.subvec(burn_in, total_time - 1);
  }

  // Step 7: Return the transposed data matrices in a list
  return Rcpp::List::create(
      Rcpp::Named("id") = id,
      Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
      Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
      Rcpp::Named("x") = x_t.t());
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-2-ou-fixed-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM2OUFixed)]]
Rcpp::List SimSSM2OUFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_l, const arma::vec& mu,
                          const arma::mat& phi, const arma::mat& sigma_l,
                          const arma::vec& nu, const arma::mat& lambda,
                          const arma::mat& theta_l, const arma::mat& gamma_y,
                          const arma::mat& gamma_eta, const Rcpp::List& x,
                          const double delta_t, const int time,
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
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();

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
    arma::vec sigma_vec = arma::vectorise(sigma_l * sigma_l.t());
    arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                        (arma::expmat(neg_phi_hashtag * delta_t) - J) *
                        sigma_vec;
    arma::mat psi_l =
        arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(num_latent_vars)) +
                 (gamma_eta * x_t.col(0));
    y.col(0) = nu + (lambda * eta.col(0)) +
               (theta_l * arma::randn(num_manifest_vars)) +
               (gamma_y * x_t.col(0));

    // Step 3.4: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars)) +
                   (gamma_eta * x_t.col(t));
      y.col(t) = nu + (lambda * eta.col(t)) +
                 (theta_l * arma::randn(num_manifest_vars)) +
                 (gamma_y * x_t.col(t));
    }

    // Step 3.5: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.6: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.7: Return the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return results
  return out;
}
// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-ou-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM2OUIVary)]]
Rcpp::List SimSSM2OUIVary(const int n, const Rcpp::List& mu0,
                          const Rcpp::List& sigma0_l, const Rcpp::List& mu,
                          const Rcpp::List& phi, const Rcpp::List& sigma_l,
                          const Rcpp::List& nu, const Rcpp::List& lambda,
                          const Rcpp::List& theta_l, const Rcpp::List& gamma_y,
                          const Rcpp::List& gamma_eta, const Rcpp::List& x,
                          const double delta_t, const int time,
                          const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  arma::vec mu0_temp = mu0[0];
  arma::vec nu_temp = nu[0];
  int num_latent_vars = mu0_temp.n_elem;
  int num_manifest_vars = nu_temp.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::mat y(num_manifest_vars, total_time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    arma::vec mu_temp = mu[i];
    arma::mat phi_temp = phi[i];
    arma::mat sigma_l_temp = sigma_l[i];
    arma::vec nu_temp = nu[i];
    arma::mat lambda_temp = lambda[i];
    arma::mat theta_l_temp = theta_l[i];
    arma::mat gamma_y_temp = gamma_y[i];
    arma::mat gamma_eta_temp = gamma_eta[i];

    // Step 3.2: Get state space parameters
    arma::mat I = arma::eye<arma::mat>(num_latent_vars, num_latent_vars);
    arma::mat J = arma::eye<arma::mat>(num_latent_vars * num_latent_vars,
                                       num_latent_vars * num_latent_vars);
    arma::mat neg_phi = -1 * phi_temp;
    // 3.2.1 beta_temp
    arma::mat beta_temp = arma::expmat(neg_phi * delta_t);  // A(Delta t)
    // 3.2.2 alpha_temp
    arma::vec alpha_temp = arma::inv(neg_phi) * (beta_temp - I) *
                           (phi_temp * mu_temp);  // b(Delta t)

    // 3.2.3 psi
    arma::mat neg_phi_hashtag = arma::kron(neg_phi, I) + arma::kron(I, neg_phi);
    arma::vec sigma_vec = arma::vectorise(sigma_l_temp * sigma_l_temp.t());
    arma::vec psi_vec = arma::inv(neg_phi_hashtag) *
                        (arma::expmat(neg_phi_hashtag * delta_t) - J) *
                        sigma_vec;
    arma::mat psi_l =
        arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(num_latent_vars)) +
                 (gamma_eta_temp * x_t.col(0));
    y.col(0) = nu_temp + (lambda_temp * eta.col(0)) +
               (theta_l_temp * arma::randn(num_manifest_vars)) +
               (gamma_y_temp * x_t.col(0));

    // Step 3.4: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha_temp + (beta_temp * eta.col(t - 1)) +
                   (psi_l * arma::randn(num_latent_vars)) +
                   (gamma_eta_temp * x_t.col(t));
      y.col(t) = nu_temp + (lambda_temp * eta.col(t)) +
                 (theta_l_temp * arma::randn(num_manifest_vars)) +
                 (gamma_y_temp * x_t.col(t));
    }

    // Step 3.5: If there is a burn-in period, remove it
    if (burn_in > 0) {
      y = y.cols(burn_in, total_time - 1);
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.6: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.7: Return the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("id") = id,
        Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_t.t());
  }

  // Step 4: Return results
  return out;
}
