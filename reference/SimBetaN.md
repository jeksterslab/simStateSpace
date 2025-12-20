# Simulate Transition Matrices from the Multivariate Normal Distribution

This function simulates random transition matrices from the multivariate
normal distribution. The function ensures that the generated transition
matrices are stationary using
[`TestStationarity()`](https://github.com/jeksterslab/simStateSpace/reference/TestStationarity.md)
with a rejection sampling approach.

## Usage

``` r
SimBetaN(
  n,
  beta,
  vcov_beta_vec_l,
  margin = 1,
  beta_lbound = NULL,
  beta_ubound = NULL,
  bound = FALSE,
  max_iter = 100000L
)
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

  Numeric scalar specifying the stationarity threshold. Values less than
  1 indicate stricter stationarity criteria.

- beta_lbound:

  Optional numeric matrix of same dim as `beta`. Use NA for no lower
  bound.

- beta_ubound:

  Optional numeric matrix of same dim as `beta`. Use NA for no upper
  bound.

- bound:

  Logical; if TRUE, resample until all elements respect bounds (NA
  bounds ignored).

- max_iter:

  Safety cap on resampling attempts per draw.

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
[`SSMMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanEta.md),
[`SSMMeanY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanY.md),
[`SimAlphaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimAlphaN.md),
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
SimBetaN(n = n, beta = beta, vcov_beta_vec_l = vcov_beta_vec_l)
#> [[1]]
#>            [,1]         [,2]        [,3]
#> [1,]  0.6851834 -0.007643585  0.03225375
#> [2,]  0.4989322  0.535273121 -0.01040462
#> [3,] -0.1017514  0.431520771  0.52115628
#> 
#> [[2]]
#>            [,1]       [,2]         [,3]
#> [1,]  0.6680390 0.02629276  0.005338959
#> [2,]  0.4718801 0.58538447 -0.013719421
#> [3,] -0.1012473 0.35777458  0.502049507
#> 
#> [[3]]
#>             [,1]       [,2]       [,3]
#> [1,]  0.72547649 0.03917744 0.01601526
#> [2,]  0.46980611 0.58514000 0.01645385
#> [3,] -0.09449514 0.39903621 0.48528219
#> 
#> [[4]]
#>             [,1]       [,2]       [,3]
#> [1,]  0.67274293 0.06733423 0.01210433
#> [2,]  0.55638168 0.60140159 0.02169502
#> [3,] -0.09678768 0.43988592 0.50541587
#> 
#> [[5]]
#>             [,1]       [,2]         [,3]
#> [1,]  0.67997152 0.05736071 -0.005371841
#> [2,]  0.53904114 0.63490134  0.009108167
#> [3,] -0.06383272 0.40650940  0.517786160
#> 
#> [[6]]
#>            [,1]         [,2]       [,3]
#> [1,]  0.6840106 -0.006206552 0.04608043
#> [2,]  0.5443914  0.569218457 0.06780036
#> [3,] -0.1077355  0.414977012 0.52529448
#> 
#> [[7]]
#>            [,1]        [,2]        [,3]
#> [1,]  0.7124042 -0.04549078  0.02242287
#> [2,]  0.4636874  0.62246809 -0.02782470
#> [3,] -0.0838195  0.36851786  0.48349206
#> 
#> [[8]]
#>            [,1]       [,2]        [,3]
#> [1,]  0.7137180 0.02645314  0.02205182
#> [2,]  0.4502689 0.62743727 -0.02979570
#> [3,] -0.1022793 0.40554066  0.46948154
#> 
#> [[9]]
#>            [,1]       [,2]        [,3]
#> [1,]  0.7166506 0.04992734 -0.02545914
#> [2,]  0.4881213 0.65434223  0.04074997
#> [3,] -0.1123248 0.37019981  0.50786389
#> 
#> [[10]]
#>             [,1]        [,2]         [,3]
#> [1,]  0.69686953 -0.02832438 -0.009376019
#> [2,]  0.44792962  0.60969799  0.026254137
#> [3,] -0.08429006  0.39618472  0.507779452
#> 
```
