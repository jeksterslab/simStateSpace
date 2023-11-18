// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-ou-fixed-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0OUFixed)]]
Rcpp::List SimSSM0OUFixed(const int n, const arma::vec& mu0, const arma::mat& sigma0_sqrt, const arma::vec& mu, const arma::mat& phi, const arma::mat& sigma_sqrt, const arma::vec& nu, const arma::mat& lambda, const arma::mat& theta_sqrt, const double delta_t, const int time, const int burn_in) {
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
    arma::mat J = arma::eye<arma::mat>(num_latent_vars * num_latent_vars, num_latent_vars * num_latent_vars);
    arma::mat neg_phi = -1 * phi;
    // 3.2.1 beta
    arma::mat beta = arma::expmat(neg_phi * delta_t);  // A(Delta t)
    // 3.2.2 alpha
    arma::vec alpha = arma::inv(neg_phi) * (beta - I) * (phi * mu);  // b(Delta t)
    
    // 3.2.3 psi
    arma::mat neg_phi_hashtag = arma::kron(neg_phi, I) + arma::kron(I, neg_phi);
    arma::vec sigma_vec = arma::vectorise(sigma_sqrt * sigma_sqrt.t());
    arma::vec psi_vec = arma::inv(neg_phi_hashtag) * (arma::expmat(neg_phi_hashtag * delta_t) - J) * sigma_vec;
    arma::mat psi_sqrt = arma::chol(arma::reshape(psi_vec, num_latent_vars, num_latent_vars));

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0 + sigma0_sqrt * arma::randn(num_latent_vars);
    y.col(0) = nu + lambda * eta.col(0) + theta_sqrt * arma::randn(num_manifest_vars);

    // Step 3.4: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + beta * eta.col(t - 1) + psi_sqrt * arma::randn(num_latent_vars);
      y.col(t) = nu + lambda * eta.col(t) + theta_sqrt * arma::randn(num_manifest_vars);
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
    out[i] = Rcpp::List::create(Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(), Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time), Rcpp::Named("id") = id);
  }

  // Step 4: Return results
  return out;
}
