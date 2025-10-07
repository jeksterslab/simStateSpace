// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-stationarity.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test Stationarity
//'
//' The function computes the eigenvalues of the input matrix `x`.
//' It checks if all eigenvalues have moduli less than 1.
//' If all eigenvalues have moduli less than 1,
//' the system is considered stationary.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric matrix.
//' @param r_target Numeric scalar specifying the stationarity threshold.
//'   Values less than 1 indicate stricter stationarity criteria.
//'
//' @examples
//' x <- matrix(
//'   data = c(0.5, 0.3, 0.2, 0.4),
//'   nrow = 2
//' )
//' TestStationarity(x)
//'
//' x <- matrix(
//'   data = c(0.9, -0.5, 0.8, 0.7),
//'   nrow = 2
//' )
//' TestStationarity(x)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test ssm
//' @export
// [[Rcpp::export]]
bool TestStationarity(const arma::mat& x, const double r_target = 1.0) {
  arma::cx_vec eigenvalues = arma::eig_gen(x);
  return arma::all(arma::abs(eigenvalues) < r_target);
}
