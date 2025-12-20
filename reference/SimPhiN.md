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
[`SimCovN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovN.md),
[`SimIotaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimIotaN.md),
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
#>            [,1]        [,2]         [,3]
#> [1,] -0.3644239 -0.00686005  0.009970126
#> [2,]  0.7239323 -0.52938142  0.003689034
#> [3,] -0.4719376  0.74081034 -0.754595776
#> 
#> [[2]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.4035570  0.03114586 -0.03271647
#> [2,]  0.8064403 -0.49974876  0.03096005
#> [3,] -0.4538086  0.71467118 -0.71241195
#> 
#> [[3]]
#>            [,1]          [,2]        [,3]
#> [1,] -0.3510828  0.0001627612 -0.06830001
#> [2,]  0.7582763 -0.4885667794 -0.01073329
#> [3,] -0.4247432  0.7258476066 -0.70662470
#> 
#> [[4]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3735633 -0.04210809 -0.08426100
#> [2,]  0.7876941 -0.47384904  0.04763003
#> [3,] -0.4729038  0.71066065 -0.69490587
#> 
#> [[5]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3434258  0.02084666 -0.004789062
#> [2,]  0.8196683 -0.55031101  0.016136026
#> [3,] -0.4922228  0.72245455 -0.674902868
#> 
#> [[6]]
#>            [,1]       [,2]        [,3]
#> [1,] -0.3461487  0.0322910  0.01262215
#> [2,]  0.7449142 -0.5253309  0.03305843
#> [3,] -0.4975340  0.7853186 -0.73998345
#> 
#> [[7]]
#>            [,1]        [,2]          [,3]
#> [1,] -0.3615692 -0.03668001  0.0005504809
#> [2,]  0.7318196 -0.49750777 -0.0274337166
#> [3,] -0.4485742  0.69805598 -0.6928824748
#> 
#> [[8]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3304681 -0.01838064 -0.04678486
#> [2,]  0.7534764 -0.50637754  0.03727028
#> [3,] -0.4956332  0.74469756 -0.68462958
#> 
#> [[9]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3975887 -0.05272651 -0.085302392
#> [2,]  0.7927197 -0.53432838 -0.004704189
#> [3,] -0.4492825  0.73642496 -0.716512771
#> 
#> [[10]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3446627  0.02023395  0.03049497
#> [2,]  0.7672337 -0.52994584 -0.01360389
#> [3,] -0.4769160  0.71338521 -0.68765222
#> 
```
