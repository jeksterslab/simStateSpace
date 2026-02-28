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
#>            [,1]       [,2]       [,3]
#> [1,]  0.6982486 0.03152077 0.02115628
#> [2,]  0.4923564 0.63225375 0.07626581
#> [3,] -0.1647269 0.38959538 0.46803899
#> 
#> [[2]]
#>            [,1]        [,2]        [,3]
#> [1,]  0.6987527 -0.04222542 0.002049507
#> [2,]  0.5262928  0.60533896 0.038227590
#> [3,] -0.1146155  0.38628058 0.525476493
#> 
#> [[3]]
#>            [,1]          [,2]        [,3]
#> [1,]  0.7055049 -0.0009637857 -0.01471781
#> [2,]  0.5391774  0.6160152581 -0.04312751
#> [3,] -0.1148600  0.4164538481  0.47274293
#> 
#> [[4]]
#>             [,1]       [,2]         [,3]
#> [1,]  0.70321232 0.03988592 5.415868e-03
#> [2,]  0.56733423 0.61210433 2.089554e-06
#> [3,] -0.09859841 0.42169502 4.799715e-01
#> 
#> [[5]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.73616728 0.006509402  0.01778616
#> [2,]  0.55736071 0.594628159 -0.01217511
#> [3,] -0.06509866 0.409108167  0.48401056
#> 
#> [[6]]
#>            [,1]       [,2]        [,3]
#> [1,]  0.6922645 0.01497701  0.02529448
#> [2,]  0.4937934 0.64608043 -0.03444810
#> [3,] -0.1307815 0.46780036  0.51240420
#> 
#> [[7]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.71618050 -0.03148214 -0.01650794
#> [2,]  0.45450922  0.62242287 -0.06306017
#> [3,] -0.07753191  0.37217530  0.51371805
#> 
#> [[8]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.69772075 0.005540657 -0.03051846
#> [2,]  0.52645314 0.622051822  0.03249946
#> [3,] -0.07256273 0.370204304  0.51665058
#> 
#> [[9]]
#>             [,1]        [,2]         [,3]
#> [1,]  0.68767516 -0.02980019  0.007863893
#> [2,]  0.54992734  0.57454086 -0.019854149
#> [3,] -0.04565777  0.44074997  0.496869531
#> 
#> [[10]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.71570994 -0.00381528 0.007779452
#> [2,]  0.47167562  0.59062398 0.005225181
#> [3,] -0.09030201  0.42625414 0.495154569
#> 
```
