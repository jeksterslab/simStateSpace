// -----------------------------------------------------------------------------
// edit .setup/cpp/simStateSpace-sim-beta-n-0.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

inline bool in_bounds_element(double x, double lb, double ub, bool has_lb,
                              bool has_ub) {
  // Treat non-finite bounds (NA/NaN/Inf) as "no bound" on that side
  bool lb_ok = true, ub_ok = true;
  if (has_lb && std::isfinite(lb)) lb_ok = (x >= lb);
  if (has_ub && std::isfinite(ub)) ub_ok = (x <= ub);
  return lb_ok && ub_ok;
}

inline bool matrix_in_bounds(const arma::mat& x, const arma::mat* lb_ptr,
                             const arma::mat* ub_ptr, bool has_lb,
                             bool has_ub) {
  if (!has_lb && !has_ub) return true;  // nothing to check
  const arma::uword nr = x.n_rows, nc = x.n_cols;
  for (arma::uword i = 0; i < nr; ++i) {
    for (arma::uword j = 0; j < nc; ++j) {
      double lb = has_lb ? (*lb_ptr)(i, j) : NA_REAL;
      double ub = has_ub ? (*ub_ptr)(i, j) : NA_REAL;
      if (!in_bounds_element(x(i, j), lb, ub, has_lb, has_ub)) return false;
    }
  }
  return true;
}
