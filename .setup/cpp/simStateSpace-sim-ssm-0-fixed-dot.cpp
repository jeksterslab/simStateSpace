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
