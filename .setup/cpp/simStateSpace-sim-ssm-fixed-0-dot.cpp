// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-fixed-0-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMFixed0)]]
Rcpp::List SimSSMFixed0(const arma::uword& n, const arma::uword& time,
                        const double delta_t, const arma::vec& mu0,
                        const arma::mat& sigma0_l, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_l,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_l) {
  // Step 1: Determine dimensions
  // int p = mu0.n_elem; // number of latent variables
  // int k = nu.n_elem; // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector
  arma::vec id_template(time, arma::fill::zeros);

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (arma::uword i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, and id variables
    arma::mat eta(mu0.n_elem, time, arma::fill::zeros);
    arma::mat y(nu.n_elem, time, arma::fill::zeros);
    arma::vec id = id_template;
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(mu0.n_elem));
    y.col(0) = nu + (lambda * eta.col(0)) + (theta_l * arma::randn(nu.n_elem));
    // Step 3.3: Data generation loop
    for (arma::uword t = 1; t < time; t++) {
      eta.col(t) =
          alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(mu0.n_elem));
      y.col(t) =
          nu + (lambda * eta.col(t)) + (theta_l * arma::randn(nu.n_elem));
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t());
  }

  // Step 4: Return results
  return output;
}
