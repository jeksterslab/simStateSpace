// -----------------------------------------------------------------------------
// edit .setup/cpp/000-forward-declarations.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List SimSSM0(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                   const arma::vec& alpha, const arma::mat& beta,
                   const arma::mat& psi_sqrt, const arma::vec& nu,
                   const arma::mat& lambda, const arma::mat& theta_sqrt,
                   const int time, const int burn_in);

Rcpp::List SimSSMVAR(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                     const arma::vec& alpha, const arma::mat& beta,
                     const arma::mat& psi_sqrt, const int time,
                     const int burn_in);

Rcpp::List SimSSMOU(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                    const arma::vec& mu, const arma::mat& phi,
                    const arma::mat& sigma_sqrt, const arma::vec& nu,
                    const arma::mat& lambda, const arma::mat& theta_sqrt,
                    const double delta_t, const int time, const int burn_in);

Rcpp::List OU2SSM(const arma::vec& mu, const arma::mat& phi,
                  const arma::mat& sigma_sqrt, const double delta_t);

Rcpp::List SimSSM0Fixed(const int n, const arma::vec& mu0,
                        const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_sqrt,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_sqrt, const int time,
                        const int burn_in);

Rcpp::List SimSSMVARFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                          const arma::mat& beta, const arma::mat& psi_sqrt,
                          const int time, const int burn_in);

Rcpp::List SimSSMOUFixed(const int n, const arma::vec& mu0,
                         const arma::mat& sigma0_sqrt, const arma::vec& mu,
                         const arma::mat& phi, const arma::mat& sigma_sqrt,
                         const arma::vec& nu, const arma::mat& lambda,
                         const arma::mat& theta_sqrt, const double delta_t,
                         const int time, const int burn_in);
