# Simulate Random Drift Matrices from the Multivariate Normal Distribution

This function simulates random drift matrices from the multivariate
normal distribution. The function ensures that the generated drift
matrices are stable using
[`TestPhi()`](https://github.com/jeksterslab/simStateSpace/reference/TestPhi.md).

## Usage

``` r
SimPhiN(
  n,
  phi,
  vcov_phi_vec_l,
  margin = 0,
  auto_ubound = 0,
  phi_lbound = NULL,
  phi_ubound = NULL,
  bound = FALSE,
  max_iter = 100000L
)
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

  Numeric scalar specifying the stability threshold for the real part of
  the eigenvalues. The default `0.0` corresponds to the imaginary axis;
  values less than `0.0` enforce a stricter stability margin.

- auto_ubound:

  Numeric scalar specifying the upper bound for the diagonal elements of
  \\\boldsymbol{\Phi}\\. Default is `0.0`, requiring all diagonal values
  to be \\\leq 0\\.

- phi_lbound:

  Optional numeric matrix of same dim as `phi`. Use NA for no lower
  bound.

- phi_ubound:

  Optional numeric matrix of same dim as `phi`. Use NA for no upper
  bound.

- bound:

  Logical; if TRUE, resample until all elements respect bounds (NA
  bounds ignored).

- max_iter:

  Safety cap on resampling attempts per draw.

## Value

Returns a list of random drift matrices.

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
phi <- matrix(
  data = c(
    -0.357, 0.771, -0.450,
    0.0, -0.511, 0.729,
    0, 0, -0.693
  ),
  nrow = 3
)
vcov_phi_vec_l <- t(chol(0.001 * diag(9)))
SimPhiN(n = n, phi = phi, vcov_phi_vec_l = vcov_phi_vec_l)
#> [[1]]
#>            [,1]         [,2]         [,3]
#> [1,] -0.3803284 -0.004704189  0.012337300
#> [2,]  0.7784250 -0.534512771 -0.003766315
#> [3,] -0.5353024  0.688890510 -0.719916012
#> 
#> [[2]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3759458 -0.01360389  0.008170644
#> [2,]  0.7553852 -0.50565222 -0.001605182
#> [3,] -0.4195050  0.73861802 -0.748566704
#> 
#> [[3]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3578406  0.01103531 -0.007005884
#> [2,]  0.7192967 -0.55382263  0.016173962
#> [3,] -0.5015185  0.75114202 -0.642153269
#> 
#> [[4]]
#>            [,1]         [,2]         [,3]
#> [1,] -0.3234106 -0.002503593  0.007110629
#> [2,]  0.8053209 -0.469841743  0.030643889
#> [3,] -0.4668201  0.737721520 -0.688656053
#> 
#> [[5]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3747053  0.01818416  0.02044528
#> [2,]  0.7173553 -0.54661089  0.02281127
#> [3,] -0.4077556  0.71871552 -0.71797204
#> 
#> [[6]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3711852 -0.07951748 -0.044538028
#> [2,]  0.7691125 -0.50222099 -0.003600132
#> [3,] -0.4155389  0.70863383 -0.636222032
#> 
#> [[7]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3352115  0.02281463 -0.02307038
#> [2,]  0.7590190 -0.57323524  0.01159043
#> [3,] -0.4828905  0.72346135 -0.71078982
#> 
#> [[8]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3548964  0.00341812  0.03168255
#> [2,]  0.7914933 -0.53605078 -0.02504568
#> [3,] -0.4606836  0.72500330 -0.75334504
#> 
#> [[9]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3221946 -0.02493372 -0.03409528
#> [2,]  0.8082950 -0.51975468 -0.01031097
#> [3,] -0.4578512  0.73499965 -0.68843429
#> 
#> [[10]]
#>            [,1]         [,2]        [,3]
#> [1,] -0.3192668 -0.008600326 -0.00254667
#> [2,]  0.7260043 -0.532294682  0.03013701
#> [3,] -0.4576226  0.731226638 -0.69282861
#> 
```
