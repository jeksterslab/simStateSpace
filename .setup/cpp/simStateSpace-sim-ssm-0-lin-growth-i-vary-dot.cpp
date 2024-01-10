// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-0-lin-growth-i-vary-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.SimSSM0LinGrowthIVary)]]
Rcpp::List SimSSM0LinGrowthIVary(const int n, const Rcpp::List& mu0,
                                 const Rcpp::List& sigma0_l,
                                 const Rcpp::List& theta_l, const int time) {
  // Step 1: Create constant vectors and matrices
  arma::mat lambda = {{1, 0}};
  arma::mat beta = {{1, 1}, {0, 1}};

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(2, time);
    arma::mat y(1, time);
    arma::vec mu0_temp = mu0[i];
    arma::mat sigma0_l_temp = sigma0_l[i];
    double theta_l_temp = theta_l[i];

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0_temp + (sigma0_l_temp * arma::randn(2));
    y.col(0) = (lambda * eta.col(0)) + (theta_l_temp * arma::randn(1));

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < time; t++) {
      eta.col(t) = (beta * eta.col(t - 1));
      y.col(t) = (lambda * eta.col(t)) + (theta_l_temp * arma::randn(1));
    }

    // Step 3.4: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.5: Save the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("y") = y.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("x") = 0, Rcpp::Named("time") = arma::regspace(0, time - 1),
        Rcpp::Named("id") = id);
  }

  // Step 4: Return the results
  return out;
}
