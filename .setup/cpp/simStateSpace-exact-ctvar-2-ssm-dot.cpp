// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-exact-ctvar-2-ssm-dot.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.CTVAR2SSM)]]
Rcpp::List CTVAR2SSM(const arma::vec& iota, const Rcpp::List& phi,
                     const arma::mat& sigma_l, const double delta_t) {
  const arma::mat a = DerivativeCompanion(phi);
  const arma::uword n = a.n_rows, q = iota.n_elem;
  arma::vec b(n, arma::fill::zeros);
  b.tail(q) = iota;
  arma::mat g(n, n, arma::fill::zeros);
  g.submat(n - q, n - q, n - 1, n - 1) = sigma_l;
  Rcpp::List out = LinSDE2SSM(b, a, g, delta_t);
  arma::mat lambda(q, n, arma::fill::zeros);
  lambda.submat(0, 0, q - 1, q - 1) = arma::eye(q, q);
  out["lambda"] = lambda;
  return out;
}
