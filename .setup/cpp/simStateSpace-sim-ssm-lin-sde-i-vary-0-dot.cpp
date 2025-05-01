#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLinSDEIVary0)]]
Rcpp::List SimSSMLinSDEIVary0(const arma::uword& n, const arma::uword& time,
                              const double delta_t, const Rcpp::List& mu0,
                              const Rcpp::List& sigma0_l,
                              const Rcpp::List& iota, const Rcpp::List& phi,
                              const Rcpp::List& sigma_l, const Rcpp::List& nu,
                              const Rcpp::List& lambda,
                              const Rcpp::List& theta_l,
                              const bool ou = false) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  arma::vec nu_i = nu[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = nu_i.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);
  arma::mat I = arma::eye(mu0_i.n_elem, mu0_i.n_elem);
  arma::mat J =
      arma::eye(mu0_i.n_elem * mu0_i.n_elem, mu0_i.n_elem * mu0_i.n_elem);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(nu_i.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec iota_i = Rcpp::as<arma::vec>(iota[i]);
    arma::mat phi_i = phi[i];
    arma::mat sigma_l_i = sigma_l[i];
    arma::vec nu_i = Rcpp::as<arma::vec>(nu[i]);
    arma::mat lambda_i = lambda[i];
    arma::mat theta_l_i = theta_l[i];

    // Step 3.3: Calculate state space parameters
    // if (ou) {
    //   iota_i = (-1 * phi_i) * iota_i;
    // }
    // arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
    // arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
    // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
    // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i; arma::mat
    // psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem,
    // mu0_i.n_elem), "lower"); arma::mat beta_i = arma::expmat(phi_i *
    // delta_t); arma::vec alpha_i = arma::inv(phi_i) * (beta_i - I) * iota_i;
    if (ou && !iota_i.is_zero()) {
      iota_i = (-1 * phi_i) * iota_i;
    }
    arma::vec alpha_i = arma::vec(mu0_i.n_elem);
    arma::mat beta_i = arma::expmat(phi_i * delta_t);
    // alpha_i = iota_i.is_zero() ? iota_i : arma::inv(phi_i) * (beta_i - I) *
    // iota_i;
    alpha_i =
        iota_i.is_zero() ? iota_i : arma::solve(phi_i, (beta_i - I) * iota_i);
    arma::mat psi_l_i = arma::mat(mu0_i.n_elem, mu0_i.n_elem, arma::fill::none);
    if (sigma_l_i.is_zero()) {
      psi_l_i = sigma_l_i;
    } else {
      arma::mat phi_hashtag_i = arma::kron(phi_i, I) + arma::kron(I, phi_i);
      arma::vec sigma_vec_i = arma::vectorise(sigma_l_i * sigma_l_i.t());
      // arma::vec psi_vec_i = arma::inv(phi_hashtag_i) *
      // (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i;
      arma::vec psi_vec_i = arma::solve(
          phi_hashtag_i,
          (arma::expmat(phi_hashtag_i * delta_t) - J) * sigma_vec_i);
      psi_l_i = arma::chol(arma::reshape(psi_vec_i, mu0_i.n_elem, mu0_i.n_elem),
                           "lower");
    }

    // Step 3.4: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem));
    y.col(0) =
        nu_i + (lambda_i * eta.col(0)) + (theta_l_i * arma::randn(nu_i.n_elem));
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem));
      y.col(t) = nu_i + (lambda_i * eta.col(t)) +
                 (theta_l_i * arma::randn(nu_i.n_elem));
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
