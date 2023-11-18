// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0)]]
Rcpp::List SimSSM0(const arma::vec& mu0, const arma::mat& sigma0_sqrt, const arma::vec& alpha, const arma::mat& beta, const arma::mat& psi_sqrt, const arma::vec& nu, const arma::mat& lambda, const arma::mat& theta_sqrt, const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;
  int num_manifest_vars = nu.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat y(num_manifest_vars, total_time);

  // Step 3: Generate initial condition
  eta.col(0) = mu0 + sigma0_sqrt * arma::randn(num_latent_vars);
  y.col(0) = nu + lambda * eta.col(0) + theta_sqrt * arma::randn(num_manifest_vars);

  // Step 4: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + beta * eta.col(t - 1) + psi_sqrt * arma::randn(num_latent_vars);
    y.col(t) = nu + lambda * eta.col(t) + theta_sqrt * arma::randn(num_manifest_vars);
  }

  // Step 5: If there is a burn-in period, remove it
  if (burn_in > 0) {
    y = y.cols(burn_in, total_time - 1);
    eta = eta.cols(burn_in, total_time - 1);
  }

  // Step 6: Return the transposed data matrices in a list
  return Rcpp::List::create(Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(), Rcpp::Named("time") = arma::regspace(1, time));
}
