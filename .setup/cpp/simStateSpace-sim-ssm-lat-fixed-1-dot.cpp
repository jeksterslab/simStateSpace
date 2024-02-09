// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-lat-fixed-1-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSMLatFixed1)]]
Rcpp::List SimSSMLatFixed1(const int n, const int time, const double delta_t,
                           const arma::vec& mu0, const arma::mat& sigma0_l,
                           const arma::vec& alpha, const arma::mat& beta,
                           const arma::mat& psi_l, const Rcpp::List& x,
                           const arma::mat& gamma) {
  // Step 1: Determine dimensions
  int p = mu0.n_elem;  // number of latent variables
  int k = p;           // number of observed variables
  arma::vec time_vec =
      arma::linspace(0, (time - 1) * delta_t, time);  // time vector

  // Step 2: Initialize the output list
  Rcpp::List output(n);

  // Step 3: Generate data per individual
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices of latent, observed, covariate, and id
    // variables
    arma::mat eta(p, time);
    arma::mat y(k, time);
    arma::mat x_i = x[i];
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);
    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + (sigma0_l * arma::randn(p)) + (gamma * x_i.col(0));
    y.col(0) = eta.col(0);
    // Step 3.3: Data generation loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = alpha + (beta * eta.col(t - 1)) + (psi_l * arma::randn(p)) +
                   (gamma * x_i.col(t));
      y.col(t) = eta.col(t);
    }
    // Step 3.4 Save results in a list
    output[i] = Rcpp::List::create(
        Rcpp::Named("id") = id, Rcpp::Named("time") = time_vec,
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = x_i.t());
  }

  // Step 4: Return results
  return output;
}
