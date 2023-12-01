// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-1-var-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM1VARIVary)]]
Rcpp::List SimSSM1VARIVary(const int n, const Rcpp::List& mu0, const Rcpp::List& sigma0_sqrt, const Rcpp::List& alpha, const Rcpp::List& beta, const Rcpp::List& psi_sqrt, const Rcpp::List& gamma_eta, const Rcpp::List& x, const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  arma::vec mu0_temp = mu0[0];
  int num_latent_vars = mu0_temp.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);
    arma::mat x_temp = x[i];
    arma::mat x_t = x_temp.t();
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_sqrt_temp = sigma0_sqrt[i];
    arma::vec alpha_temp = alpha[i];
    arma::mat beta_temp = beta[i];
    arma::mat psi_sqrt_temp = psi_sqrt[i];
    arma::mat gamma_eta_temp = gamma_eta[i];

    // Step 3.2: Generate initial condition    
    eta.col(0) = mu0_temp + sigma0_sqrt_temp * arma::randn(num_latent_vars) + gamma_eta_temp * x_t.col(0);

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha_temp + beta_temp * eta.col(t - 1) + psi_sqrt_temp * arma::randn(num_latent_vars) + gamma_eta_temp * x_t.col(t);
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      eta = eta.cols(burn_in, total_time - 1);
      x_t = x_t.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(Rcpp::Named("y") = eta.t(), Rcpp::Named("eta") = eta.t(), Rcpp::Named("x") = x_t.t(), Rcpp::Named("time") = arma::regspace(0, time - 1), Rcpp::Named("id") = id);
  }

  // Step 4: Return the results
  return out;
}
