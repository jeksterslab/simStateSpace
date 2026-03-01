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
#>             [,1]        [,2]        [,3]
#> [1,]  0.63527312 -0.01040462 -0.03196101
#> [2,]  0.53152077  0.62115628 -0.02811988
#> [3,] -0.06774625  0.47626581  0.49875269
#> 
#> [[2]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.68538447 -0.01371942  0.02547649
#> [2,]  0.45777458  0.60204951 -0.03019389
#> [3,] -0.09466104  0.43822759  0.50550486
#> 
#> [[3]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.68514000 0.01645385 -0.02725707
#> [2,]  0.49903621 0.58528219  0.05638168
#> [3,] -0.08398474 0.35687249  0.50321232
#> 
#> [[4]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.70140159 0.02169502 -0.02002848
#> [2,]  0.53988592 0.60541587  0.03904114
#> [3,] -0.08789567 0.40000209  0.53616728
#> 
#> [[5]]
#>            [,1]        [,2]        [,3]
#> [1,]  0.7349013 0.009108167 -0.01598944
#> [2,]  0.5065094 0.617786160  0.04439135
#> [3,] -0.1053718 0.387824885  0.49226451
#> 
#> [[6]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.66921846 0.06780036  0.01240420
#> [2,]  0.51497701 0.62529448 -0.03631255
#> [3,] -0.05391957 0.36555190  0.51618050
#> 
#> [[7]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.72246809 -0.0278247  0.01371805
#> [2,]  0.46851786  0.5834921 -0.04973111
#> [3,] -0.07757713  0.3369398  0.49772075
#> 
#> [[8]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.72743727 -0.0297957  0.01665058
#> [2,]  0.50554066  0.5694815 -0.01187867
#> [3,] -0.07794818  0.4324995  0.48767516
#> 
#> [[9]]
#>            [,1]       [,2]         [,3]
#> [1,]  0.7543422 0.04074997 -0.003130469
#> [2,]  0.4701998 0.60786389 -0.052070380
#> [3,] -0.1254591 0.38014585  0.515709935
#> 
#> [[10]]
#>            [,1]       [,2]         [,3]
#> [1,]  0.7096980 0.02625414 -0.004845431
#> [2,]  0.4961847 0.60777945  0.014382456
#> [3,] -0.1093760 0.40522518  0.479894529
#> 
```
