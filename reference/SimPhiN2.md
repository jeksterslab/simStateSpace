# Simulate Random Drift Matrices from the Multivariate Normal Distribution and Project to Hurwitz

This function simulates random dirft matrices from the multivariate
normal distribution then projects each draw to the Hurwitz-stable region
using
[`ProjectToHurwitz()`](https://github.com/jeksterslab/simStateSpace/reference/ProjectToHurwitz.md).

## Usage

``` r
SimPhiN2(n, phi, vcov_phi_vec_l, margin = 0.001)
```

## Arguments

- n:

  Positive integer. Number of replications.

- phi:

  Numeric matrix. The drift matrix (\\\boldsymbol{\Phi}\\).

- vcov_phi_vec_l:

  Numeric matrix. Cholesky factorization (`t(chol(vcov_phi_vec))`) of
  the sampling variance-covariance matrix of \\\mathrm{vec} \left(
  \boldsymbol{\Phi} \right)\\.

- margin:

  Positive numeric. Target buffer so that the spectral abscissa is \\\le
  -\text{margin}\\ (default `1e-3`).

## Value

Returns a list of random drift matrices.

## See also

Other Simulation of State Space Models Data Functions:
[`LinSDE2SSM()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDE2SSM.md),
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
[`SimMVN()`](https://github.com/jeksterslab/simStateSpace/reference/SimMVN.md),
[`SimMuN()`](https://github.com/jeksterslab/simStateSpace/reference/SimMuN.md),
[`SimNuN()`](https://github.com/jeksterslab/simStateSpace/reference/SimNuN.md),
[`SimPhiN()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN.md),
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
n <- 10
phi <- matrix(
  data = c(
    -0.357, 0.771, -0.450,
    0.0, -0.511, 0.729,
    0, 0, -0.693
  ),
  nrow = 3
)
vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
SimPhiN2(n = n, phi = phi, vcov_phi_vec_l = vcov_phi_vec_l)
#> [[1]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3374115 -0.03631547 -0.03522833
#> [2,]  0.7829520 -0.51424701 -0.02122889
#> [3,] -0.4593276  0.71458045 -0.72769393
#> 
#> [[2]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3639378 -0.02178019  0.03261244
#> [2,]  0.7968141 -0.51485211  0.00289716
#> [3,] -0.4635005  0.69815458 -0.66238104
#> 
#> [[3]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3942297  0.02094066 -0.02158031
#> [2,]  0.7984825 -0.53096047  0.03250162
#> [3,] -0.4819357  0.74745080 -0.70935891
#> 
#> [[4]]
#>            [,1]       [,2]        [,3]
#> [1,] -0.3657105  0.0333104 -0.02528062
#> [2,]  0.8121938 -0.5232222 -0.04990546
#> [3,] -0.4275957  0.7362539 -0.64011685
#> 
#> [[5]]
#>            [,1]         [,2]        [,3]
#> [1,] -0.3276867  0.002224122  0.02757431
#> [2,]  0.7527464 -0.564012342 -0.02370732
#> [3,] -0.4262468  0.724154016 -0.70016594
#> 
#> [[6]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3434785 -0.03717281 -0.026153021
#> [2,]  0.8416886 -0.56042143  0.001184558
#> [3,] -0.4619536  0.76184134 -0.731455246
#> 
#> [[7]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3765914 -0.01058848 -0.04617385
#> [2,]  0.8166669 -0.54811423 -0.01558777
#> [3,] -0.4811700  0.68307304 -0.66308147
#> 
#> [[8]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.2827693 -0.04749673 -0.02039162
#> [2,]  0.7541351 -0.53188626  0.08655416
#> [3,] -0.4355784  0.71807782 -0.64235630
#> 
#> [[9]]
#>            [,1]         [,2]        [,3]
#> [1,] -0.3355311 -0.009941172 -0.02135158
#> [2,]  0.7472695 -0.540460856 -0.05944341
#> [3,] -0.4138522  0.734597964 -0.71644021
#> 
#> [[10]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.4417425 -0.01756875  0.001208828
#> [2,]  0.7840805 -0.48789034  0.017578063
#> [3,] -0.4818337  0.74712907 -0.708353500
#> 
```
