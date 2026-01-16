# Convert Parameters from the Linear Stochastic Differential Equation Model to State Space Model Parameterization

This function converts parameters from the linear stochastic
differential equation model to state space model parameterization.

## Usage

``` r
LinSDE2SSM(iota, phi, sigma_l, delta_t)
```

## Arguments

- iota:

  Numeric vector. An unobserved term that is constant over time
  (\\\boldsymbol{\iota}\\).

- phi:

  Numeric matrix. The drift matrix which represents the rate of change
  of the solution in the absence of any random fluctuations
  (\\\boldsymbol{\Phi}\\).

- sigma_l:

  Numeric matrix. Cholesky factorization (`t(chol(sigma))`) of the
  covariance matrix of volatility or randomness in the process
  (\\\boldsymbol{\Sigma}\\).

- delta_t:

  Numeric. Time interval (\\\Delta_t\\).

## Value

Returns a list of state space parameters:

- `alpha`: Numeric vector. Vector of constant values for the dynamic
  model (\\\boldsymbol{\alpha}\\).

- `beta`: Numeric matrix. Transition matrix relating the values of the
  latent variables from the previous time point to the current time
  point. (\\\boldsymbol{\beta}\\).

- `psi_l`: Numeric matrix. Cholesky factorization (`t(chol(psi))`) of
  the process noise covariance matrix \\\boldsymbol{\Psi}\\.

## Details

Let the linear stochastic equation model be given by \$\$ \mathrm{d}
\boldsymbol{\eta}\_{i, t} = \left( \boldsymbol{\iota} +
\boldsymbol{\Phi} \boldsymbol{\eta}\_{i, t} \right) \mathrm{d} t +
\boldsymbol{\Sigma}^{\frac{1}{2}} \mathrm{d} \mathbf{W}\_{i, t} \$\$ for
individual \\i\\ and time \\t\\. The discrete-time state space model
given below represents the discrete-time solution for the linear
stochastic differential equation. \$\$ \boldsymbol{\eta}\_{i,
t\_{{l\_{i}}}} = \boldsymbol{\alpha}\_{\Delta t\_{{l\_{i}}}} +
\boldsymbol{\beta}\_{\Delta t\_{{l\_{i}}}} \boldsymbol{\eta}\_{i,
t\_{l\_{i} - 1}} + \boldsymbol{\zeta}\_{i, t\_{{l\_{i}}}}, \quad
\mathrm{with} \quad \boldsymbol{\zeta}\_{i, t\_{{l\_{i}}}} \sim
\mathcal{N} \left( \mathbf{0}, \boldsymbol{\Psi}\_{\Delta t\_{{l\_{i}}}}
\right) \$\$ with \$\$ \boldsymbol{\beta}\_{\Delta t\_{{l\_{i}}}} =
\exp{ \left( \Delta t \boldsymbol{\Phi} \right) }, \$\$

\$\$ \boldsymbol{\alpha}\_{\Delta t\_{{l\_{i}}}} =
\boldsymbol{\Phi}^{-1} \left( \boldsymbol{\beta} - \mathbf{I}\_{p}
\right) \boldsymbol{\iota}, \quad \mathrm{and} \$\$

\$\$ \mathrm{vec} \left( \boldsymbol{\Psi}\_{\Delta t\_{{l\_{i}}}}
\right) = \left\[ \left( \boldsymbol{\Phi} \otimes \mathbf{I}\_{p}
\right) + \left( \mathbf{I}\_{p} \otimes \boldsymbol{\Phi} \right)
\right\] \left\[ \exp \left( \left\[ \left( \boldsymbol{\Phi} \otimes
\mathbf{I}\_{p} \right) + \left( \mathbf{I}\_{p} \otimes
\boldsymbol{\Phi} \right) \right\] \Delta t \right) - \mathbf{I}\_{p
\times p} \right\] \mathrm{vec} \left( \boldsymbol{\Sigma} \right) \$\$
where \\t\\ denotes continuous-time processes that can be defined by any
arbitrary time point, \\t\_{l\_{i}}\\ the \\l^\mathrm{th}\\ observed
measurement occassion for individual \\i\\, \\p\\ the number of latent
variables and \\\Delta t\\ the time interval.

## References

Harvey, A. C. (1990). Forecasting, structural time series models and the
Kalman filter. Cambridge University Press.
[doi:10.1017/cbo9781107049994](https://doi.org/10.1017/cbo9781107049994)

## See also

Other Simulation of State Space Models Data Functions:
[`LinSDECovEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDECovEta.md),
[`LinSDECovY()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDECovY.md),
[`LinSDEMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEMeanEta.md),
[`LinSDEMeanY()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEMeanY.md),
[`ProjectToHurwitz()`](https://github.com/jeksterslab/simStateSpace/reference/ProjectToHurwitz.md),
[`ProjectToStability()`](https://github.com/jeksterslab/simStateSpace/reference/ProjectToStability.md),
[`SSMCovEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMCovEta.md),
[`SSMCovY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMCovY.md),
[`SSMInterceptEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMInterceptEta.md),
[`SSMInterceptY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMInterceptY.md),
[`SSMMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanEta.md),
[`SSMMeanY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanY.md),
[`SimAlphaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimAlphaN.md),
[`SimBetaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN.md),
[`SimBetaN2()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN2.md),
[`SimBetaNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaNCovariate.md),
[`SimCovDiagN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovDiagN.md),
[`SimCovN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovN.md),
[`SimIotaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimIotaN.md),
[`SimNuN()`](https://github.com/jeksterslab/simStateSpace/reference/SimNuN.md),
[`SimPhiN()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN.md),
[`SimPhiN2()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN2.md),
[`SimPhiNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiNCovariate.md),
[`SimSSMFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMFixed.md),
[`SimSSMIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMIVary.md),
[`SimSSMLinGrowth()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinGrowth.md),
[`SimSSMLinGrowthIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinGrowthIVary.md),
[`SimSSMLinSDEFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEFixed.md),
[`SimSSMLinSDEIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEIVary.md),
[`SimSSMOUFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUFixed.md),
[`SimSSMOUIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUIVary.md),
[`SimSSMVARFixed()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMVARFixed.md),
[`SimSSMVARIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMVARIVary.md),
[`SpectralRadius()`](https://github.com/jeksterslab/simStateSpace/reference/SpectralRadius.md),
[`TestPhi()`](https://github.com/jeksterslab/simStateSpace/reference/TestPhi.md),
[`TestPhiHurwitz()`](https://github.com/jeksterslab/simStateSpace/reference/TestPhiHurwitz.md),
[`TestStability()`](https://github.com/jeksterslab/simStateSpace/reference/TestStability.md),
[`TestStationarity()`](https://github.com/jeksterslab/simStateSpace/reference/TestStationarity.md)

## Author

Ivan Jacob Agaloos Pesigan

## Examples

``` r
p <- 2
iota <- c(0.317, 0.230)
phi <- matrix(
  data = c(
   -0.10,
   0.05,
   0.05,
   -0.10
 ),
 nrow = p
)
sigma <- matrix(
  data = c(
    2.79,
    0.06,
    0.06,
    3.27
  ),
  nrow = p
)
sigma_l <- t(chol(sigma))
delta_t <- 0.10

LinSDE2SSM(
  iota = iota,
  phi = phi,
  sigma_l = sigma_l,
  delta_t = delta_t
)
#> $alpha
#>            [,1]
#> [1,] 0.03159928
#> [2,] 0.02296420
#> 
#> $beta
#>            [,1]       [,2]
#> [1,] 0.99006221 0.00495027
#> [2,] 0.00495027 0.99006221
#> 
#> $psi_l
#>            [,1]      [,2]
#> [1,] 0.52560735 0.0000000
#> [2,] 0.01414641 0.5688463
#> 
```
