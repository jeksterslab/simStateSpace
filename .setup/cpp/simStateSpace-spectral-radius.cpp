// -----------------------------------------------------------------------------
// .setup/cpp/simStateSpace-spectral-radius.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Spectral Radius
//'
//' Computes the spectral radius of a square matrix,
//' defined as the maximum modulus (absolute value) of its eigenvalues.
//' The spectral radius is often used to assess the stability of systems
//' such as vector autoregressive (VAR) models:
//' a system is considered stationary
//' if the spectral radius of its transition matrix is strictly less than 1.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric square matrix.
//'
//' @return Numeric value representing the spectral radius of `x`.
//'
//' @examples
//' # Matrix with eigenvalues less than 1
//' x <- matrix(
//'   data = c(
//'     0.5, 0.3,
//'     0.2, 0.4
//'   ),
//'   nrow = 2
//' )
//' SpectralRadius(x)
//'
//' # Matrix with eigenvalues greater than 1
//' y <- matrix(
//'   data = c(
//'     1.2, 0.3,
//'     0.4, 0.9
//'   ),
//'   nrow = 2
//' )
//' SpectralRadius(y)
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace stability ssm
//' @export
// [[Rcpp::export]]
double SpectralRadius(const arma::mat& x) {
  arma::cx_vec ev = arma::eig_gen(x);
  return arma::abs(ev).max();
}
