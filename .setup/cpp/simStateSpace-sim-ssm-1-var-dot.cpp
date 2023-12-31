// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-var-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1VAR)]]
Rcpp::List SimSSM1VAR(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                      const arma::vec& alpha, const arma::mat& beta,
                      const arma::mat& psi_sqrt, const arma::mat& gamma_eta,
                      const arma::mat& x, const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;

  // Step 2: Create matrices to store simulated data
  arma::mat eta(num_latent_vars, total_time);
  arma::mat x_t = x.t();
  arma::vec id(total_time, arma::fill::ones);

  // Step 3: Generate initial condition
  eta.col(0) =
      mu0 + sigma0_sqrt * arma::randn(num_latent_vars) + gamma_eta * x_t.col(0);

  // Step 4: Simulate state space model data using a loop
  for (int t = 1; t < total_time; t++) {
    eta.col(t) = alpha + beta * eta.col(t - 1) +
                 psi_sqrt * arma::randn(num_latent_vars) +
                 gamma_eta * x_t.col(t);
  }

  // Step 5: If there is a burn-in period, remove it
  if (burn_in > 0) {
    eta = eta.cols(burn_in, total_time - 1);
    x_t = x_t.cols(burn_in, total_time - 1);
    id = id.subvec(burn_in, total_time - 1);
  }

  // Step 6: Return the transposed data matrices in a list
  return Rcpp::List::create(Rcpp::Named("y") = eta.t(),
                            Rcpp::Named("eta") = eta.t(),
                            Rcpp::Named("x") = x_t.t(),
                            Rcpp::Named("time") = arma::regspace(0, time - 1),
                            Rcpp::Named("id") = id);
}
