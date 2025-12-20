# Simulate Covariance Matrices from the Multivariate Normal Distribution

This function simulates random covariance matrices from the multivariate
normal distribution. The function ensures that the generated covariance
matrices are positive semi-definite.

## Usage

``` r
SimCovN(n, sigma, vcov_sigma_vech_l)
```

## Arguments

- n:

  Positive integer. Number of replications.

- sigma:

  Numeric matrix. The covariance matrix (\\\boldsymbol{\Sigma}\\).

- vcov_sigma_vech_l:

  Numeric matrix. Cholesky factorization (`t(chol(vcov_sigma_vech))`) of
  the sampling variance-covariance matrix of \\\mathrm{vech} \left(
  \boldsymbol{\Sigma} \right)\\.

## Value

Returns a list of random covariance matrices.

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
[`SSMMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanEta.md),
[`SSMMeanY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanY.md),
[`SimAlphaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimAlphaN.md),
[`SimBetaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN.md),
[`SimBetaN2()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN2.md),
[`SimBetaNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaNCovariate.md),
[`SimCovDiagN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovDiagN.md),
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
n <- 10
sigma <- matrix(
  data = c(
    1.0, 0.5, 0.3,
    0.5, 1.0, 0.4,
    0.3, 0.4, 1.0
  ),
  nrow = 3
)
vcov_sigma_vech_l <- t(
  chol(
    0.001 * diag(3 * (3 + 1) / 2)
  )
)
SimCovN(
  n = n,
  sigma = sigma,
  vcov_sigma_vech_l = vcov_sigma_vech_l
)
#> [[1]]
#>           [,1]      [,2]      [,3]
#> [1,] 1.0009730 0.5796971 0.2974445
#> [2,] 0.5796971 1.0318089 0.4512592
#> [3,] 0.2974445 0.4512592 0.9680744
#> 
#> [[2]]
#>           [,1]      [,2]      [,3]
#> [1,] 1.0099154 0.4502194 0.2802876
#> [2,] 0.4502194 1.0110675 0.3887407
#> [3,] 0.2802876 0.3887407 0.9849274
#> 
#> [[3]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9747801 0.4793178 0.3656572
#> [2,] 0.4793178 0.9663087 0.3904455
#> [3,] 0.3656572 0.3904455 0.9793431
#> 
#> [[4]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9594177 0.5274351 0.3060267
#> [2,] 0.5274351 1.0358792 0.3990601
#> [3,] 0.3060267 0.3990601 0.9455439
#> 
#> [[5]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9357604 0.4534859 0.3150958
#> [2,] 0.4534859 1.0011291 0.3767855
#> [3,] 0.3150958 0.3767855 1.0697661
#> 
#> [[6]]
#>           [,1]      [,2]      [,3]
#> [1,] 1.0185430 0.5475251 0.4089177
#> [2,] 0.5475251 0.9857585 0.3988875
#> [3,] 0.4089177 0.3988875 1.0114199
#> 
#> [[7]]
#>           [,1]      [,2]      [,3]
#> [1,] 1.0212965 0.4333275 0.3001881
#> [2,] 0.4333275 0.9745739 0.4073115
#> [3,] 0.3001881 0.4073115 1.0023112
#> 
#> [[8]]
#>           [,1]      [,2]      [,3]
#> [1,] 1.0106630 0.5007951 0.3483639
#> [2,] 0.5007951 1.0119281 0.4389563
#> [3,] 0.3483639 0.4389563 1.0194549
#> 
#> [[9]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9588815 0.4522809 0.3379359
#> [2,] 0.4522809 0.9952336 0.4135109
#> [3,] 0.3379359 0.4135109 1.0350359
#> 
#> [[10]]
#>           [,1]      [,2]      [,3]
#> [1,] 1.0295486 0.5430591 0.3211406
#> [2,] 0.5430591 0.9875579 0.4086113
#> [3,] 0.3211406 0.4086113 1.0350278
#> 
```
