# Simulate Transition Matrices from the Multivariate Normal Distribution and Project to Stability

This function simulates random transition matrices from the multivariate
normal distribution then projects each draw to the stability region
using
[`ProjectToStability()`](https://github.com/jeksterslab/simStateSpace/reference/ProjectToStability.md).

## Usage

``` r
SimBetaN2(n, beta, vcov_beta_vec_l, margin = 0.98, tol = 1e-12)
```

## Arguments

- n:

  Positive integer. Number of replications.

- beta:

  Numeric matrix. The transition matrix (\\\boldsymbol{\beta}\\).

- vcov_beta_vec_l:

  Numeric matrix. Cholesky factorization (`t(chol(vcov_beta_vec))`) of
  the sampling variance-covariance matrix of \\\mathrm{vec} \left(
  \boldsymbol{\beta} \right)\\.

- margin:

  Double in \\(0, 1)\\. Target upper bound for the spectral radius
  (default = 0.98).

- tol:

  Small positive double added to the denominator in the scaling factor
  to avoid division by zero (default = 1e-12).

## Value

Returns a list of random transition matrices.

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
n <- 10
beta <- matrix(
  data = c(
    0.7, 0.5, -0.1,
    0.0, 0.6, 0.4,
    0, 0, 0.5
  ),
  nrow = 3
)
vcov_beta_vec_l <- t(chol(0.001 * diag(9)))
SimBetaN2(n = n, beta = beta, vcov_beta_vec_l = vcov_beta_vec_l)
#> [[1]]
#>            [,1]       [,2]       [,3]
#> [1,]  0.6798945 0.02458508 0.01605138
#> [2,]  0.4796956 0.58802190 0.04630413
#> [3,] -0.1223069 0.42123165 0.51984513
#> 
#> [[2]]
#>            [,1]        [,2]       [,3]
#> [1,]  0.6923103 0.007366145 0.01225458
#> [2,]  0.4663594 0.634151590 0.02122083
#> [3,] -0.1197942 0.452563278 0.55052020
#> 
#> [[3]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.68444080 -0.04965351  0.04652801
#> [2,]  0.50400738  0.63639433 -0.03453232
#> [3,] -0.04750596  0.40179583  0.51801095
#> 
#> [[4]]
#>             [,1]        [,2]       [,3]
#> [1,]  0.72412758 -0.04984393 0.04302717
#> [2,]  0.53740576  0.60905180 0.02610376
#> [3,] -0.09181022  0.31579241 0.49660088
#> 
#> [[5]]
#>             [,1]        [,2]       [,3]
#> [1,]  0.69688925 -0.04475027 0.05322493
#> [2,]  0.49305056  0.62310965 0.07649833
#> [3,] -0.08791394  0.41979430 0.47440904
#> 
#> [[6]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.67124464 0.01205901 0.005107774
#> [2,]  0.51673866 0.56243996 0.050158536
#> [3,] -0.09526582 0.45001538 0.434261085
#> 
#> [[7]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.70827280 0.01117277  0.01925151
#> [2,]  0.51146601 0.65195525 -0.03432050
#> [3,] -0.06432483 0.37508497  0.50595960
#> 
#> [[8]]
#>             [,1]         [,2]         [,3]
#> [1,]  0.68137923 -0.008413454 0.0004621614
#> [2,]  0.48806017  0.630368234 0.0142915462
#> [3,] -0.06537565  0.426496814 0.4857642965
#> 
#> [[9]]
#>            [,1]        [,2]        [,3]
#> [1,]  0.7146258 -0.03574804 -0.02434587
#> [2,]  0.4831519  0.61346426  0.01219682
#> [3,] -0.1133236  0.43964171  0.53603814
#> 
#> [[10]]
#>             [,1]      [,2]        [,3]
#> [1,]  0.67388074 0.0144374 -0.03827979
#> [2,]  0.54816789 0.5689039  0.00764590
#> [3,] -0.09329593 0.4728070  0.52755878
#> 
```
