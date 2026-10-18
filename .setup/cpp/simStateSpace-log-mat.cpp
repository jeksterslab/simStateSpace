// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-log-mat.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Matrix Logarithm
//'
//' Computes the real-valued principal matrix logarithm of a square numeric
//' matrix using the Armadillo linear algebra library.
//'
//' @param x A square numeric matrix.
//' @param tol A positive numeric scalar giving the tolerance used to determine
//'   whether the imaginary component of the matrix logarithm is negligible.
//'   The default is `1e-10`.
//'
//' @return
//' A numeric matrix containing the real part of the matrix logarithm of `x`.
//'
//' @examples
//' a <- matrix(
//'   data = c(
//'     -0.40,  0.10,
//'     -0.05, -0.30
//'   ),
//'   nrow = 2,
//'   ncol = 2
//' )
//'
//' b <- ExpMat(x = a)
//' b
//'
//' a <- LogMat(x = b)
//' a
//'
//' @family Simulation of State Space Models Data Functions
//' @keywords simStateSpace linsde
//' @export
// [[Rcpp::export]]
arma::mat LogMat(const arma::mat& x, const double tol = 1e-10) {
  if (x.n_rows != x.n_cols) {
    Rcpp::stop("`x` must be a square matrix.");
  }
  if (!x.is_finite()) {
    Rcpp::stop("`x` must contain only finite values.");
  }
  if (!std::isfinite(tol)) {
    Rcpp::stop("`tol` must be finite.");
  }
  if (tol <= 0.0) {
    Rcpp::stop("`tol` must be positive.");
  }
  arma::cx_mat out;
  const bool success = arma::logmat(out, x);
  if (!success) {
    Rcpp::stop("Armadillo failed to compute the matrix logarithm.");  // # nocov
  }
  const double imag_norm = arma::norm(arma::imag(out), "inf");
  if (!std::isfinite(imag_norm) || imag_norm > tol) {
    Rcpp::stop("Matrix logarithm has a non-negligible imaginary component.");
  }
  return arma::real(out);
}
