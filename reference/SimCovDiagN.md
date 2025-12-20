# Simulate Diagonal Covariance Matrices from the Multivariate Normal Distribution

This function simulates random diagonal covariance matrices from the
multivariate normal distribution. The function ensures that the
generated covariance matrices are positive semi-definite.

## Usage

``` r
SimCovDiagN(n, sigma_diag, vcov_sigma_diag_l)
```

## Arguments

- n:

  Positive integer. Number of replications.

- sigma_diag:

  Numeric matrix. The covariance matrix (\\\boldsymbol{\Sigma}\\).

- vcov_sigma_diag_l:

  Numeric matrix. Cholesky factorization (`t(chol(vcov_sigma_vech))`) of
  the sampling variance-covariance matrix of \\\mathrm{vech} \left(
  \boldsymbol{\Sigma} \right)\\.

## Value

Returns a list of random diagonal covariance matrices.

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
n <- 10
sigma_diag <- c(1, 1, 1)
vcov_sigma_diag_l <- t(chol(0.001 * diag(3)))
SimCovDiagN(
  n = n,
  sigma_diag = sigma_diag,
  vcov_sigma_diag_l = vcov_sigma_diag_l
)
#> [[1]]
#>          [,1]     [,2]     [,3]
#> [1,] 1.018961 0.000000 0.000000
#> [2,] 0.000000 1.034873 0.000000
#> [3,] 0.000000 0.000000 1.024973
#> 
#> [[2]]
#>           [,1]     [,2]      [,3]
#> [1,] 0.9286775 0.000000 0.0000000
#> [2,] 0.0000000 1.012595 0.0000000
#> [3,] 0.0000000 0.000000 0.9977968
#> 
#> [[3]]
#>          [,1]      [,2]     [,3]
#> [1,] 1.028824 0.0000000 0.000000
#> [2,] 0.000000 0.9096426 0.000000
#> [3,] 0.000000 0.0000000 1.021352
#> 
#> [[4]]
#>           [,1]      [,2]     [,3]
#> [1,] 0.9742215 0.0000000 0.000000
#> [2,] 0.0000000 0.9844067 0.000000
#> [3,] 0.0000000 0.0000000 1.028971
#> 
#> [[5]]
#>          [,1]      [,2]      [,3]
#> [1,] 1.023703 0.0000000 0.0000000
#> [2,] 0.000000 0.9572151 0.0000000
#> [3,] 0.000000 0.0000000 0.9597157
#> 
#> [[6]]
#>           [,1]      [,2]     [,3]
#> [1,] 0.9954098 0.0000000 0.000000
#> [2,] 0.0000000 0.9975186 0.000000
#> [3,] 0.0000000 0.0000000 1.023891
#> 
#> [[7]]
#>           [,1]    [,2]     [,3]
#> [1,] 0.9898989 0.00000 0.000000
#> [2,] 0.0000000 1.00651 0.000000
#> [3,] 0.0000000 0.00000 1.043984
#> 
#> [[8]]
#>          [,1]      [,2]      [,3]
#> [1,] 1.007064 0.0000000 0.0000000
#> [2,] 0.000000 0.9440822 0.0000000
#> [3,] 0.000000 0.0000000 0.9471224
#> 
#> [[9]]
#>           [,1]     [,2]    [,3]
#> [1,] 0.9638273 0.000000 0.00000
#> [2,] 0.0000000 1.017432 0.00000
#> [3,] 0.0000000 0.000000 1.02906
#> 
#> [[10]]
#>           [,1]      [,2]    [,3]
#> [1,] 0.9846425 0.0000000 0.00000
#> [2,] 0.0000000 0.9959778 0.00000
#> [3,] 0.0000000 0.0000000 1.02866
#> 
```
