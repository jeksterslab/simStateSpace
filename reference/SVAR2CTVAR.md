# Convert SVAR Parameters to CTVAR Parameters

This function converts structural vector autoregressive (SVAR)
parameters to continuous-time vector autoregressive (CTVAR) parameters
using the SVAR to CTVAR transformation in Table 1, block C, of Chow et
al. (2023).

## Usage

``` r
SVAR2CTVAR(alpha_star, beta_0_star, beta_star, psi_star_l, delta_t)
```

## Arguments

- alpha_star:

  Numeric vector. SVAR intercept vector \\\boldsymbol{\alpha}^{\ast}\\.
  Use a vector of zeroes when the SVAR has no intercept.

- beta_0_star:

  Numeric matrix. Contemporaneous SVAR matrix
  \\\boldsymbol{\beta}\_0^{\ast}\\.

- beta_star:

  List. A list of numeric matrices containing the SVAR lag matrices
  \\\boldsymbol{\beta}\_1^{\ast}, \ldots,
  \boldsymbol{\beta}\_p^{\ast}\\.

- psi_star_l:

  Numeric matrix. Lower Cholesky factorization (`t(chol(psi_star))`) of
  the SVAR innovation covariance matrix \\\boldsymbol{\Psi}^{\ast}\\.

- delta_t:

  Numeric. Time interval.

## Value

Returns a named list with `iota`, `phi`, `sigma_l`, and `s`.

## Details

Let \\p\\ be the SVAR order. The SVAR to CTVAR transformation solves for
\\S_0, \ldots, S\_{p - 1}\\ using \$\$ \boldsymbol{\beta}\_0^{\ast} =
I + \sum\_{m = 0}^{p} S_m \Delta t^{-m} \$\$ and \$\$
\boldsymbol{\beta}\_j^{\ast} = (-1)^j \sum\_{m = j}^{p} {m \choose j}
S_m \Delta t^{-m}, \quad j = 1, \ldots, p, \$\$ where \\S_p = -I\\.

The returned `phi` is the first-order companion-form drift matrix for
the CTVAR process. The returned `iota` contains `alpha_star` in the
final block. The returned `sigma_l` contains `psi_star_l` in the final
block.

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
[`VAR2CTVAR()`](https://github.com/jeksterslab/simStateSpace/reference/VAR2CTVAR.md),
[`VAR2SVAR()`](https://github.com/jeksterslab/simStateSpace/reference/VAR2SVAR.md)

## Author

Ivan Jacob Agaloos Pesigan

## Examples

``` r
# Example 1: SVAR(1) to CTVAR(1)
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

svar <- VAR2SVAR(
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
  delta_t = 1.0
)

out <- SVAR2CTVAR(
  alpha_star = svar$alpha_star,
  beta_0_star = svar$beta_0_star,
  beta_star = svar$beta_star,
  psi_star_l = svar$psi_star_l,
  delta_t = 1.0
)

out$iota
#>            [,1]
#> [1,]  0.2291022
#> [2,] -0.1238390
out$phi
#>            [,1]       [,2]
#> [1,] -0.1145511  0.0619195
#> [2,]  0.0619195 -0.1145511
out$sigma_l
#>           [,1]      [,2]
#> [1,] 1.1034883 0.0000000
#> [2,] 0.1132674 0.9783552

# Example 2: SVAR(2) to CTVAR(2)
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

svar <- VAR2SVAR(
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
  delta_t = 1.0
)

out <- SVAR2CTVAR(
  alpha_star = svar$alpha_star,
  beta_0_star = svar$beta_0_star,
  beta_star = svar$beta_star,
  psi_star_l = svar$psi_star_l,
  delta_t = 1.0
)

out$iota
#>           [,1]
#> [1,] 0.0000000
#> [2,] 0.0000000
#> [3,] 0.5789474
#> [4,] 0.5263158
out$phi
#>            [,1]       [,2]       [,3]       [,4]
#> [1,]  0.0000000  0.0000000  1.0000000  0.0000000
#> [2,]  0.0000000  0.0000000  0.0000000  1.0000000
#> [3,] -3.5263158 -0.3684211 -0.7368421 -0.6842105
#> [4,] -0.1754386 -5.1228070 -1.5789474 -0.8947368
out$sigma_l
#>      [,1] [,2]     [,3]     [,4]
#> [1,]    0    0 0.000000 0.000000
#> [2,]    0    0 0.000000 0.000000
#> [3,]    0    0 4.600373 0.000000
#> [4,]    0    0 3.211430 4.883756

# For a bivariate CTVAR(2), the companion-form state
# has dimension 4 by 4.
dim(out$phi)
#> [1] 4 4

# Example 3: SVAR(2) with zero intercept and zero innovation covariance
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

svar <- VAR2SVAR(
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
  delta_t = 1.0
)

out <- SVAR2CTVAR(
  alpha_star = svar$alpha_star,
  beta_0_star = svar$beta_0_star,
  beta_star = svar$beta_star,
  psi_star_l = svar$psi_star_l,
  delta_t = 1.0
)

out$iota
#>      [,1]
#> [1,]    0
#> [2,]    0
#> [3,]    0
#> [4,]    0
out$phi
#>            [,1]        [,2]       [,3]       [,4]
#> [1,]  0.0000000  0.00000000  1.0000000  0.0000000
#> [2,]  0.0000000  0.00000000  0.0000000  1.0000000
#> [3,] -3.4092742 -0.06048387  0.3770161 -0.1411290
#> [4,] -0.2318548 -4.51209677 -0.5745968  0.4717742
out$sigma_l
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    0    0    0    0
#> [3,]    0    0    0    0
#> [4,]    0    0    0    0
```
