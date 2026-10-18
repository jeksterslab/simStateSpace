// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-companion-mat-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.CompanionMat)]]
arma::mat CompanionMat(const Rcpp::List& x) {
  const R_xlen_t order_r = x.size();
  Rcpp::NumericMatrix x_1_r(x[0]);
  const arma::uword q = static_cast<arma::uword>(x_1_r.nrow());
  const arma::uword order = static_cast<arma::uword>(order_r);
  const arma::uword n = q * order;
  arma::mat out = arma::zeros(n, n);
  const arma::mat identity = arma::eye(q, q);
  for (arma::uword i = 0; i < order; i++) {
    Rcpp::NumericMatrix x_i_r(x[i]);
    const arma::mat x_i(x_i_r.begin(), x_i_r.nrow(), x_i_r.ncol(), false);
    out.submat(0, i * q, q - 1, ((i + 1) * q) - 1) = x_i;
  }
  for (arma::uword i = 1; i < order; i++) {
    out.submat(i * q, (i - 1) * q, ((i + 1) * q) - 1, (i * q) - 1) = identity;
  }
  return out;
}
