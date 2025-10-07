// -----------------------------------------------------------------------------
// edit .setup/cpp/000-forward-declarations.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

bool TestPhi(const arma::mat& phi, const double a_target,
             const double auto_ubound);

bool TestStability(const arma::mat& x, const double a_target);

bool TestStationarity(const arma::mat& x, const double r_target);

double SpectralRadius(const arma::mat& x);

arma::mat ProjectToStability(const arma::mat& x, const double margin,
                             const double tol);

double SpectralAbscissa(const arma::mat& x);

arma::mat ProjectToHurwitz(const arma::mat& x, const double margin);
