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