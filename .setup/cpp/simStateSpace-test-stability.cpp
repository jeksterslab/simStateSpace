// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-test-stability.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Test Stability
//'
//' The function computes the eigenvalues of the input matrix `x`.
//' It checks if the real part of all eigenvalues is negative.
//' If all eigenvalues have negative real parts,
//' the system is considered stable.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric matrix.
//' @param margin Numeric scalar specifying the stability threshold
//'   for the real part of the eigenvalues.
//'   The default `0.0` corresponds to the imaginary axis;
//'   values less than `0.0` enforce a stricter stability margin.
//'
//' @examples
//' x <- matrix(
//'   data = c(
//'     -0.357, 0.771, -0.450,
//'     0.0, -0.511, 0.729,
//'     0, 0, -0.693
//'   ),
//'   nrow = 3
//' )
//' TestStability(x)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace test linsde
//' @export
// [[Rcpp::export]]
bool TestStability(const arma::mat& x, const double margin = 0.0) {
  arma::cx_vec eigenvalues = arma::eig_gen(x);
  return arma::all(arma::real(eigenvalues) < margin);
}
