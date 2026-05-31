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
#>            [,1]        [,2]        [,3]
#> [1,] -0.3617891 -0.01156665 -0.04753395
#> [2,]  0.7871360 -0.50014870  0.03229100
#> [3,] -0.4319029  0.70291418 -0.70733089
#> 
#> [[2]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3443779  0.04254781  0.001425817
#> [2,]  0.8040584 -0.51556917 -0.036680008
#> [3,] -0.4969834  0.68981961 -0.679507768
#> 
#> [[3]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3564495 -0.01791724 -0.04563322
#> [2,]  0.7435663 -0.48446806 -0.01838064
#> [3,] -0.4498825  0.71147641 -0.68837754
#> 
#> [[4]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.4037849 -0.00914248  0.000717545
#> [2,]  0.8082703 -0.55158870 -0.052726506
#> [3,] -0.4416296  0.75071973 -0.716328379
#> 
#> [[5]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.4423024 -0.04010949 -0.02691601
#> [2,]  0.7662958 -0.49866270  0.02023395
#> [3,] -0.4735128  0.72523368 -0.71194584
#> 
#> [[6]]
#>            [,1]         [,2]         [,3]
#> [1,] -0.3265050  0.009618024 -0.055566704
#> [2,]  0.7573961 -0.502829356 -0.002825644
#> [3,] -0.4446522  0.727394818 -0.693840568
#> 
#> [[7]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.4085185  0.02214202  0.05084673
#> [2,]  0.7820353 -0.51800588  0.01491232
#> [3,] -0.4928226  0.74517396 -0.65941060
#> 
#> [[8]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3738201  0.00872152  0.004343947
#> [2,]  0.7684964 -0.50388937  0.008644179
#> [3,] -0.4088417  0.75964389 -0.710705326
#> 
#> [[9]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3147556 -0.01028448 -0.02497204
#> [2,]  0.7891842 -0.49055472  0.02514259
#> [3,] -0.4856109  0.75181127 -0.70718519
#> 
#> [[10]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3225389 -0.02036617  0.05677797
#> [2,]  0.6914825 -0.55553803  0.01008559
#> [3,] -0.4412210  0.72539987 -0.67121153
#> 
```
