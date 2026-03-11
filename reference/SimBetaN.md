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
#>            [,1]        [,2]         [,3]
#> [1,]  0.7113137 -0.06108355 -0.001098922
#> [2,]  0.4973354  0.55895700  0.048238135
#> [3,] -0.1128323  0.40398333  0.471195115
#> 
#> [[2]]
#>             [,1]        [,2]        [,3]
#> [1,]  0.68236604 -0.01742727  0.01244699
#> [2,]  0.49408692  0.62516893 -0.03988741
#> [3,] -0.07494344  0.34096968  0.52048583
#> 
#> [[3]]
#>            [,1]        [,2]         [,3]
#> [1,]  0.7111349 -0.06607776 -0.001751393
#> [2,]  0.4447285  0.58518340 -0.007643585
#> [3,] -0.1274092  0.39893216  0.435273121
#> 
#> [[4]]
#>             [,1]       [,2]         [,3]
#> [1,]  0.73225375 0.07626581 -0.001247313
#> [2,]  0.48959538 0.56803899  0.026292756
#> [3,] -0.07884372 0.37188012  0.485384475
#> 
#> [[5]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.70533896 0.03822759 0.005504865
#> [2,]  0.48628058 0.62547649 0.039177440
#> [3,] -0.09795049 0.36980611 0.485140003
#> 
#> [[6]]
#>            [,1]        [,2]        [,3]
#> [1,]  0.7160153 -0.04312751 0.003212323
#> [2,]  0.5164538  0.57274293 0.067334232
#> [3,] -0.1147178  0.45638168 0.501401591
#> 
#> [[7]]
#>             [,1]         [,2]       [,3]
#> [1,]  0.71210433 2.089554e-06 0.03616728
#> [2,]  0.52169502 5.799715e-01 0.05736071
#> [3,] -0.09458413 4.390411e-01 0.53490134
#> 
#> [[8]]
#>             [,1]        [,2]         [,3]
#> [1,]  0.69462816 -0.01217511 -0.007735489
#> [2,]  0.50910817  0.58401056 -0.006206552
#> [3,] -0.08221384  0.44439135  0.469218457
#> 
#> [[9]]
#>             [,1]       [,2]        [,3]
#> [1,]  0.74608043 -0.0344481  0.01618050
#> [2,]  0.56780036  0.6124042 -0.04549078
#> [3,] -0.07470552  0.3636874  0.52246809
#> 
#> [[10]]
#>            [,1]        [,2]         [,3]
#> [1,]  0.7224229 -0.06306017 -0.002279254
#> [2,]  0.4721753  0.61371805  0.026453144
#> [3,] -0.1165079  0.35026889  0.527437274
#> 
```
