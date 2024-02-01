// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lin-sde-i-vary-2-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLinSDEIVary2)]]
Rcpp::List SimSSMLinSDEIVary2(
    const int n, const int time, const double delta_t, const Rcpp::List& mu0,
    const Rcpp::List& sigma0_l, const Rcpp::List& gamma, const Rcpp::List& phi,
    const Rcpp::List& sigma_l, const Rcpp::List& nu, const Rcpp::List& lambda,
    const Rcpp::List& theta_l, const Rcpp::List& x, const Rcpp::List& gamma_eta,
    const Rcpp::List& gamma_y, const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  int p = mu0_i.n_elem;  // number of latent variables
  int k = nu_i.n_elem;   // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::mat I = arma::eye<arma::mat>(p, p);
  arma::mat J = arma::eye<arma::mat>(p * p, p * p);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = mu0[i];
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec gamma_i = gamma[i];
    arma::mat phi_i = phi[i];
    arma::mat sigma_l_i = sigma_l[i];
    arma::vec nu_i = nu[i];
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];
    arma::mat gamma_eta_i = gamma_eta[i];
    arma::mat gamma_y_i = gamma_y[i];

    // Step 3.3: Calculate state space parameters
    if (ou) {
      gamma_i = phi_i * gamma_i;
    }
    arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
                          (arma::expmat(phi_hashtag_i * delta_t) - J) *
                          sigma_vec_i;
    arma::mat psi_l_i = arma::chol(arma::reshape(psi_vec_i, p, p), "lower");
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * gamma_i;

    // Step 3.4: Generate initial condition
    eta.col(0) =
        mu0_i + (sigma0_l_i * arma::randn(p)) + (gamma_eta_i * x_i.col(0));
    y.col(0) = nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(k)) +
               (gamma_y_i * x_i.col(0));
    // Step 3.4: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(p)) + (gamma_eta_i * x_i.col(t));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) + (theta_l_i * arma::randn(k)) +
                 (gamma_y_i * x_i.col(t));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
