// -----------------------------------------------------------------------------
// .setup/cpp/simStateSpace-spectral-abscissa.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Spectral Abscissa
//'
//' Returns the maximum real part of the eigenvalues of a square matrix.
//' For continuous-time stability (Hurwitz),
//' a matrix is stable if the spectral abscissa is strictly less than 0.
//'
//' @author Ivan Jacob Agaloos Pesigan
//'
//' @param x Numeric square matrix.
//'
//' @return Numeric value \eqn{\alpha(x) = \max \Re(\lambda_i(x))}.
//'
//' @examples
//' # Hurwitz-stable (spectral abscissa < 0):
//' x <- matrix(
//'   data = c(
//'     -0.5, -0.2,
//'      1.0, -0.3
//'   ),
//'   nrow = 2
//' )
//' SpectralAbscissa(x = x) # < 0
//'
//' # Unstable (spectral abscissa > 0):
//' x <- matrix(
//'   data = c(
//'      0.10, 0.50,
//'     -0.40, 0.20
//'   ),
//'   nrow = 2
//' )
//' SpectralAbscissa(x = x) # > 0
//'
//' @keywords simStateSpace stability linsde
//' @export
// [[Rcpp::export]]
double SpectralAbscissa(const arma::mat& x) {
  arma::cx_vec ev = arma::eig_gen(x);
  return arma::max(arma::real(ev));
}
