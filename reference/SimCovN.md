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
#> [1,] 0.9742215 0.4844067 0.3289707
#> [2,] 0.4844067 1.0065213 0.4237028
#> [3,] 0.3289707 0.4237028 0.9572151
#> 
#> [[2]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9597157 0.4643785 0.2954098
#> [2,] 0.4643785 0.9975186 0.4238913
#> [3,] 0.2954098 0.4238913 0.9484577
#> 
#> [[3]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9898989 0.5065099 0.3439838
#> [2,] 0.5065099 0.9841791 0.4070641
#> [3,] 0.3439838 0.4070641 0.9440822
#> 
#> [[4]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9471224 0.5097071 0.2638273
#> [2,] 0.5097071 1.0174320 0.4290599
#> [3,] 0.2638273 0.4290599 0.9734041
#> 
#> [[5]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9846425 0.4959778 0.3286602
#> [2,] 0.4959778 0.9919400 0.3803007
#> [3,] 0.3286602 0.3803007 1.0068224
#> 
#> [[6]]
#>           [,1]      [,2]      [,3]
#> [1,] 1.0009730 0.5796971 0.2974445
#> [2,] 0.5796971 1.0318089 0.4512592
#> [3,] 0.2974445 0.4512592 0.9680744
#> 
#> [[7]]
#>           [,1]      [,2]      [,3]
#> [1,] 1.0099154 0.4502194 0.2802876
#> [2,] 0.4502194 1.0110675 0.3887407
#> [3,] 0.2802876 0.3887407 0.9849274
#> 
#> [[8]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9747801 0.4793178 0.3656572
#> [2,] 0.4793178 0.9663087 0.3904455
#> [3,] 0.3656572 0.3904455 0.9793431
#> 
#> [[9]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9594177 0.5274351 0.3060267
#> [2,] 0.5274351 1.0358792 0.3990601
#> [3,] 0.3060267 0.3990601 0.9455439
#> 
#> [[10]]
#>           [,1]      [,2]      [,3]
#> [1,] 0.9357604 0.4534859 0.3150958
#> [2,] 0.4534859 1.0011291 0.3767855
#> [3,] 0.3150958 0.3767855 1.0697661
#> 
```
