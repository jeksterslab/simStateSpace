// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-ssm-var-fixed.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Simulate Data from a Vector Autoregressive Model
//' using a State Space Model Parameterization
//' for n > 1 Individuals (Fixed Parameters)
//'
//' This function simulates data from a vector autoregressive model
//' using a state space model parameterization
//' for `n > 1` individuals.
//' In this model,
//' the parameters are invariant across individuals.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @inheritParams SimSSM0Fixed
//' @inherit SimSSM0Fixed return
//' @inherit SimSSM0 references
//'
//' @examples
//' # prepare parameters
//' set.seed(42)
//' k <- 3
//' iden <- diag(k)
//' iden_sqrt <- chol(iden)
//' null_vec <- rep(x = 0, times = k)
//' n <- 5
//' mu0 <- null_vec
//' sigma0_sqrt <- iden_sqrt
//' alpha <- null_vec
//' beta <- diag(x = 0.5, nrow = k)
//' psi_sqrt <- iden_sqrt
//' time <- 50
//' burn_in <- 0
//'
//' ssm <- SimSSMVARFixed(
//'   n = n,
//'   mu0 = mu0,
//'   sigma0_sqrt = sigma0_sqrt,
//'   alpha = alpha,
//'   beta = beta,
//'   psi_sqrt = psi_sqrt,
//'   time = time,
//'   burn_in = burn_in
//' )
//'
//' str(ssm)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace sim
//' @export
// [[Rcpp::export]]
Rcpp::List SimSSMVARFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                          const arma::mat& beta, const arma::mat& psi_sqrt,
                          const int time, const int burn_in) {
  // Step 1: Determine indices
  int total_time = time + burn_in;
  int num_latent_vars = mu0.n_elem;

  // Step 2: Create a list of length n
  Rcpp::List out(n);

  // Step 3: Loop to iterate over n
  for (int i = 0; i < n; i++) {
    // Step 3.1: Create matrices to store simulated data
    arma::mat eta(num_latent_vars, total_time);

    // Step 3.2: Generate initial condition
    eta.col(0) = mu0 + sigma0_sqrt * arma::randn(num_latent_vars);

    // Step 3.3: Simulate state space model data using a loop
    for (int t = 1; t < total_time; t++) {
      eta.col(t) = alpha + beta * eta.col(t - 1) +
                   psi_sqrt * arma::randn(num_latent_vars);
    }

    // Step 3.4: If there is a burn-in period, remove it
    if (burn_in > 0) {
      eta = eta.cols(burn_in, total_time - 1);
    }

    // Step 3.5: Create a vector of ID numbers of length time
    arma::vec id(time, arma::fill::zeros);
    id.fill(i + 1);

    // Step 3.6: Return the transposed data matrices in a list
    out[i] = Rcpp::List::create(
        Rcpp::Named("y") = eta.t(), Rcpp::Named("eta") = eta.t(),
        Rcpp::Named("time") = arma::regspace(1, time), Rcpp::Named("id") = id);
  }

  // Step 4: Return the results
  return out;
}
