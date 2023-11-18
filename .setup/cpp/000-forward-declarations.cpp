// -----------------------------------------------------------------------------
// edit .setup/cpp/000-forward-declarations.cpp
// Ivan Jacob Agaloos Pesigan
// -----------------------------------------------------------------------------

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
Rcpp::List OU2SSM(const arma::vec& mu, const arma::mat& phi,
                  const arma::mat& sigma_sqrt, const double delta_t);

Rcpp::List SimSSM0(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                   const arma::vec& alpha, const arma::mat& beta,
                   const arma::mat& psi_sqrt, const arma::vec& nu,
                   const arma::mat& lambda, const arma::mat& theta_sqrt,
                   const int time, const int burn_in);

Rcpp::List SimSSM0Fixed(const int n, const arma::vec& mu0,
                        const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_sqrt,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_sqrt, const int time,
                        const int burn_in);

Rcpp::List SimSSM0Vary(const int n, const Rcpp::List& mu0,
                       const Rcpp::List& sigma0_sqrt, const Rcpp::List& alpha,
                       const Rcpp::List& beta, const Rcpp::List& psi_sqrt,
                       const Rcpp::List& nu, const Rcpp::List& lambda,
                       const Rcpp::List& theta_sqrt, const int time,
                       const int burn_in);

Rcpp::List SimSSM0VAR(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                      const arma::vec& alpha, const arma::mat& beta,
                      const arma::mat& psi_sqrt, const int time,
                      const int burn_in);

Rcpp::List SimSSM0VARFixed(const int n, const arma::vec& mu0,
                           const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                           const arma::mat& beta, const arma::mat& psi_sqrt,
                           const int time, const int burn_in);

Rcpp::List SimSSM0VARVary(const int n, const Rcpp::List& mu0,
                          const Rcpp::List& sigma0_sqrt,
                          const Rcpp::List& alpha, const Rcpp::List& beta,
                          const Rcpp::List& psi_sqrt, const int time,
                          const int burn_in);

Rcpp::List SimSSM0OU(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                     const arma::vec& mu, const arma::mat& phi,
                     const arma::mat& sigma_sqrt, const arma::vec& nu,
                     const arma::mat& lambda, const arma::mat& theta_sqrt,
                     const double delta_t, const int time, const int burn_in);

Rcpp::List SimSSM0OUFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_sqrt, const arma::vec& mu,
                          const arma::mat& phi, const arma::mat& sigma_sqrt,
                          const arma::vec& nu, const arma::mat& lambda,
                          const arma::mat& theta_sqrt, const double delta_t,
                          const int time, const int burn_in);

Rcpp::List SimSSM0OUVary(const int n, const Rcpp::List& mu0,
                         const Rcpp::List& sigma0_sqrt, const Rcpp::List& mu,
                         const Rcpp::List& phi, const Rcpp::List& sigma_sqrt,
                         const Rcpp::List& nu, const Rcpp::List& lambda,
                         const Rcpp::List& theta_sqrt, const double delta_t,
                         const int time, const int burn_in);

Rcpp::List SimSSM1(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                   const arma::vec& alpha, const arma::mat& beta,
                   const arma::mat& psi_sqrt, const arma::vec& nu,
                   const arma::mat& lambda, const arma::mat& theta_sqrt,
                   const arma::mat& gamma_eta, const arma::mat& x,
                   const int time, const int burn_in);

Rcpp::List SimSSM1Fixed(const int n, const arma::vec& mu0,
                        const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_sqrt,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_sqrt, const arma::mat& gamma_eta,
                        const Rcpp::List& x, const int time, const int burn_in);

Rcpp::List SimSSM1Vary(const int n, const Rcpp::List& mu0,
                       const Rcpp::List& sigma0_sqrt, const Rcpp::List& alpha,
                       const Rcpp::List& beta, const Rcpp::List& psi_sqrt,
                       const Rcpp::List& nu, const Rcpp::List& lambda,
                       const Rcpp::List& theta_sqrt,
                       const Rcpp::List& gamma_eta, const Rcpp::List& x,
                       const int time, const int burn_in);

Rcpp::List SimSSM1VAR(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                      const arma::vec& alpha, const arma::mat& beta,
                      const arma::mat& psi_sqrt, const arma::mat& gamma_eta,
                      const arma::mat& x, const int time, const int burn_in);

Rcpp::List SimSSM1VARFixed(const int n, const arma::vec& mu0,
                           const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                           const arma::mat& beta, const arma::mat& psi_sqrt,
                           arma::mat& gamma_eta, const Rcpp::List& x,
                           const int time, const int burn_in);

Rcpp::List SimSSM1VARVary(const int n, const Rcpp::List& mu0,
                          const Rcpp::List& sigma0_sqrt,
                          const Rcpp::List& alpha, const Rcpp::List& beta,
                          const Rcpp::List& psi_sqrt,
                          const Rcpp::List& gamma_eta, const Rcpp::List& x,
                          const int time, const int burn_in);

Rcpp::List SimSSM1OU(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                     const arma::vec& mu, const arma::mat& phi,
                     const arma::mat& sigma_sqrt, const arma::vec& nu,
                     const arma::mat& lambda, const arma::mat& theta_sqrt,
                     const arma::mat& gamma_eta, const arma::mat& x,
                     const double delta_t, const int time, const int burn_in);

Rcpp::List SimSSM1OUFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_sqrt, const arma::vec& mu,
                          const arma::mat& phi, const arma::mat& sigma_sqrt,
                          const arma::vec& nu, const arma::mat& lambda,
                          const arma::mat& theta_sqrt,
                          const arma::mat& gamma_eta, const Rcpp::List& x,
                          const double delta_t, const int time,
                          const int burn_in);

Rcpp::List SimSSM1OUVary(const int n, const Rcpp::List& mu0,
                         const Rcpp::List& sigma0_sqrt, const Rcpp::List& mu,
                         const Rcpp::List& phi, const Rcpp::List& sigma_sqrt,
                         const Rcpp::List& nu, const Rcpp::List& lambda,
                         const Rcpp::List& theta_sqrt,
                         const Rcpp::List& gamma_eta, const Rcpp::List& x,
                         const double delta_t, const int time,
                         const int burn_in);

Rcpp::List SimSSM2(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                   const arma::vec& alpha, const arma::mat& beta,
                   const arma::mat& psi_sqrt, const arma::vec& nu,
                   const arma::mat& lambda, const arma::mat& theta_sqrt,
                   const arma::mat& gamma_y, const arma::mat& gamma_eta,
                   const arma::mat& x, const int time, const int burn_in);

Rcpp::List SimSSM2Fixed(const int n, const arma::vec& mu0,
                        const arma::mat& sigma0_sqrt, const arma::vec& alpha,
                        const arma::mat& beta, const arma::mat& psi_sqrt,
                        const arma::vec& nu, const arma::mat& lambda,
                        const arma::mat& theta_sqrt, const arma::mat& gamma_y,
                        const arma::mat& gamma_eta, const Rcpp::List& x,
                        const int time, const int burn_in);

Rcpp::List SimSSM2Vary(const int n, const Rcpp::List& mu0,
                       const Rcpp::List& sigma0_sqrt, const Rcpp::List& alpha,
                       const Rcpp::List& beta, const Rcpp::List& psi_sqrt,
                       const Rcpp::List& nu, const Rcpp::List& lambda,
                       const Rcpp::List& theta_sqrt, const Rcpp::List& gamma_y,
                       const Rcpp::List& gamma_eta, const Rcpp::List& x,
                       const int time, const int burn_in);

Rcpp::List SimSSM2OU(const arma::vec& mu0, const arma::mat& sigma0_sqrt,
                     const arma::vec& mu, const arma::mat& phi,
                     const arma::mat& sigma_sqrt, const arma::vec& nu,
                     const arma::mat& lambda, const arma::mat& theta_sqrt,
                     const arma::mat& gamma_y, const arma::mat& gamma_eta,
                     const arma::mat& x, const double delta_t, const int time,
                     const int burn_in);

Rcpp::List SimSSM2OUFixed(const int n, const arma::vec& mu0,
                          const arma::mat& sigma0_sqrt, const arma::vec& mu,
                          const arma::mat& phi, const arma::mat& sigma_sqrt,
                          const arma::vec& nu, const arma::mat& lambda,
                          const arma::mat& theta_sqrt, const arma::mat& gamma_y,
                          const arma::mat& gamma_eta, const Rcpp::List& x,
                          const double delta_t, const int time,
                          const int burn_in);

Rcpp::List SimSSM2OUVary(const int n, const Rcpp::List& mu0,
                         const Rcpp::List& sigma0_sqrt, const Rcpp::List& mu,
                         const Rcpp::List& phi, const Rcpp::List& sigma_sqrt,
                         const Rcpp::List& nu, const Rcpp::List& lambda,
                         const Rcpp::List& theta_sqrt,
                         const Rcpp::List& gamma_y, const Rcpp::List& gamma_eta,
                         const Rcpp::List& x, const double delta_t,
                         const int time, const int burn_in);
