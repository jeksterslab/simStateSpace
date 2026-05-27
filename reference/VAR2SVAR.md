# Convert VAR Parameters to SVAR Parameters

This function converts vector autoregressive (VAR) parameters to
structural vector autoregressive (SVAR) parameters using the VAR to SVAR
transformation in Table 1, block D, of Chow et al. (2023).

## Usage

``` r
VAR2SVAR(alpha, beta, psi_l, delta_t)
```

## Arguments

- alpha:

  Numeric vector. Vector of constant values for the dynamic model
  (\\\boldsymbol{\alpha}\\).

- beta:

  List. A list of numeric matrices containing the VAR transition
  matrices \\\boldsymbol{\beta}\_1, \ldots, \boldsymbol{\beta}\_p\\.

- psi_l:

  Numeric matrix. Cholesky factorization (`t(chol(psi))`) of the
  covariance matrix of the process noise (\\\boldsymbol{\Psi}\\).

- delta_t:

  Numeric. Time interval.

## Value

Returns a named list with `alpha_star`, `beta_0_star`, `beta_star`, and
`psi_star_l`.

## Details

Let \\p\\ be the VAR order and let \\\boldsymbol{\beta}\_p\\ be the
final VAR lag matrix. The VAR to SVAR transformation uses \$\$ M = I -
\boldsymbol{\beta}\_0^{\ast} = (-1)^{p + 1} \Delta t^{-p}
\boldsymbol{\beta}\_p^{-1} \$\$ where \\\boldsymbol{\beta}\_0^{\ast}\\
is the contemporaneous SVAR matrix.

The SVAR lag matrices are \$\$ \boldsymbol{\beta}\_j^{\ast} =
M\boldsymbol{\beta}\_j, \quad j = 1, \ldots, p. \$\$

The SVAR innovation covariance matrix is \$\$ \boldsymbol{\Psi}^{\ast} =
M\boldsymbol{\Psi}M'. \$\$

The SVAR intercept is \$\$ \boldsymbol{\alpha}^{\ast} =
M\boldsymbol{\alpha}. \$\$

## References

Chow, S.-M., Losardo, D., Park, J. J., & Molenaar, P. C. M. (2023).
Continuous-time dynamic models: Connections to structural equation
models and other discrete-time models. In R. H. Hoyle (Ed.), *Handbook
of structural equation modeling* (2nd ed.). The Guilford Press.

## See also

Other Simulation of State Space Models Data Functions:
[`LinSDE2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDE2SSM.md),
[`LinSDECovEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDECovEta.md),
[`LinSDECovY()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDECovY.md),
[`LinSDEInterceptEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEInterceptEta.md),
[`LinSDEInterceptY()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEInterceptY.md),
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
[`SVAR2CTVAR()`](https://github.com/jeksterslab/simStateSpace/reference/SVAR2CTVAR.md),
[`SimAlphaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimAlphaN.md),
[`SimBetaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN.md),
[`SimBetaN2()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN2.md),
[`SimBetaNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaNCovariate.md),
[`SimCovDiagN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovDiagN.md),
[`SimCovN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovN.md),
[`SimIotaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimIotaN.md),
[`SimMVN()`](https://github.com/jeksterslab/simStateSpace/reference/SimMVN.md),
[`SimMuN()`](https://github.com/jeksterslab/simStateSpace/reference/SimMuN.md),
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
[`TestStationarity()`](https://github.com/jeksterslab/simStateSpace/reference/TestStationarity.md),
[`VAR2CTVAR()`](https://github.com/jeksterslab/simStateSpace/reference/VAR2CTVAR.md)

## Author

Ivan Jacob Agaloos Pesigan

## Examples

``` r
# Example 1: VAR(1)
alpha <- c(0.20, -0.10)

beta <- list(
  matrix(
    data = c(
      0.90, 0.05,
      0.05, 0.90
    ),
    nrow = 2
  )
)

psi <- matrix(
  data = c(
    1.00, 0.20,
    0.20, 0.80
  ),
  nrow = 2
)

psi_l <- t(chol(psi))

out <- VAR2SVAR(
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
  delta_t = 1.0
)

out$alpha_star
#>            [,1]
#> [1,]  0.2291022
#> [2,] -0.1238390
out$beta_0_star
#>            [,1]       [,2]
#> [1,] -0.1145511  0.0619195
#> [2,]  0.0619195 -0.1145511
out$beta_star[[1]]
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
out$psi_star_l
#>           [,1]      [,2]
#> [1,] 1.1034883 0.0000000
#> [2,] 0.1132674 0.9783552

# For VAR(1) and delta_t = 1,
# beta_star[[1]] should be close to the identity matrix.
out$beta_star[[1]] - diag(2)
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0

# Example 2: VAR(2)
alpha <- c(0.10, 0.05)

beta <- list(
  matrix(
    data = c(
      0.50, 0.10,
      0.05, 0.40
    ),
    nrow = 2
  ),
  matrix(
    data = c(
      -0.20, 0.05,
       0.03, -0.15
    ),
    nrow = 2
  )
)

psi <- matrix(
  data = c(
    0.70, 0.10,
    0.10, 0.60
  ),
  nrow = 2
)

psi_l <- t(chol(psi))

out <- VAR2SVAR(
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
  delta_t = 1.0
)

out$alpha_star
#>           [,1]
#> [1,] 0.5789474
#> [2,] 0.5263158
out$beta_0_star
#>           [,1]      [,2]
#> [1,] -4.263158 -1.052632
#> [2,] -1.754386 -6.017544
out$beta_star[[1]]
#>          [,1]      [,2]
#> [1,] 2.736842 0.6842105
#> [2,] 1.578947 2.8947368
out$beta_star[[2]]
#>      [,1]          [,2]
#> [1,]   -1 -2.775558e-17
#> [2,]    0 -1.000000e+00
out$psi_star_l
#>          [,1]     [,2]
#> [1,] 4.600373 0.000000
#> [2,] 3.211430 4.883756

# Example 3: VAR(2) with zero intercept and zero innovation covariance
alpha <- c(0.00, 0.00)

beta <- list(
  matrix(
    data = c(
      0.40, 0.05,
      0.02, 0.30
    ),
    nrow = 2
  ),
  matrix(
    data = c(
      -0.25, 0.04,
       0.01, -0.20
    ),
    nrow = 2
  )
)

psi_l <- matrix(
  data = 0.0,
  nrow = 2,
  ncol = 2
)

out <- VAR2SVAR(
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
  delta_t = 1.0
)

out$alpha_star
#>      [,1]
#> [1,]    0
#> [2,]    0
out$beta_0_star
#>            [,1]       [,2]
#> [1,] -3.0322581 -0.2016129
#> [2,] -0.8064516 -4.0403226
out$beta_star[[1]]
#>           [,1]     [,2]
#> [1,] 1.6229839 0.141129
#> [2,] 0.5745968 1.528226
out$beta_star[[2]]
#>      [,1] [,2]
#> [1,]   -1    0
#> [2,]    0   -1
out$psi_star_l
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
```
