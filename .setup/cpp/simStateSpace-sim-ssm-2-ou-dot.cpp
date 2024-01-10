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
      Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
      Rcpp::Named("x") = x_t.t(),
      Rcpp::Named("time") = arma::linspace(0, (time - 1) * delta_t, time),
      Rcpp::Named("id") = id);
}
