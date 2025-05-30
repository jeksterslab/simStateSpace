// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-i-vary-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatIVary0)]]
Rcpp::List SimSSMLatIVary0(const arma::uword& n, const arma::uword& time,
                           const double delta_t, const Rcpp::List& mu0,
                           const Rcpp::List& sigma0_l, const Rcpp::List& alpha,
                           const Rcpp::List& beta, const Rcpp::List& psi_l) {
  // Step 1: Determine dimensions
  arma::vec mu0_i = mu0[0];
  // int p = mu0_i.n_elem; // number of latent variables
  // int k = p; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0_i.n_elem, time, arma::fill::zeros);
    arma::mat y(mu0_i.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);

    // Step 3.2: Extract the ith parameter
    arma::vec mu0_i = Rcpp::as<arma::vec>(mu0[i]);
    arma::mat sigma0_l_i = sigma0_l[i];
    arma::vec alpha_i = Rcpp::as<arma::vec>(alpha[i]);
    arma::mat beta_i = beta[i];
    arma::mat psi_l_i = psi_l[i];

    // Step 3.3: Generate initial condition
    eta.col(0) = mu0_i + (sigma0_l_i * arma::randn(mu0_i.n_elem));
    y.col(0) = eta.col(0);
    // Step 3.4: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) = alpha_i + (beta_i * eta.col(t - 1)) +
                   (psi_l_i * arma::randn(mu0_i.n_elem));
      y.col(t) = eta.col(t);
    }
    // Step 3.5 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
