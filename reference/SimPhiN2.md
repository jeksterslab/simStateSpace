# Simulate Random Drift Matrices from the Multivariate Normal Distribution and Project to Hurwitz

This function simulates random dirft matrices from the multivariate
normal distribution then projects each draw to the Hurwitz-stable region
using
[`ProjectToHurwitz()`](https://github.com/jeksterslab/simStateSpace/reference/ProjectToHurwitz.md).

## Usage

``` r
SimPhiN2(n, phi, vcov_phi_vec_l, margin = 0.001)
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

  Positive numeric. Target buffer so that the spectral abscissa is \\\le
  -\text{margin}\\ (default `1e-3`).

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
[`SimPhiN()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN.md),
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
SimPhiN2(n = n, phi = phi, vcov_phi_vec_l = vcov_phi_vec_l)
#> [[1]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.4125667 -0.05170327 -0.04282263
#> [2,]  0.7681744 -0.56251846  0.02214202
#> [3,] -0.4508406  0.74003531 -0.70000588
#> 
#> [[2]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3061533  0.03432088  0.04115826
#> [2,]  0.7859123 -0.52782012  0.00872152
#> [3,] -0.4164106  0.72649641 -0.68588937
#> 
#> [[3]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3526561 -0.05364475 -0.03561089
#> [2,]  0.7796442 -0.46875557 -0.01028448
#> [3,] -0.4677053  0.74718416 -0.67255472
#> 
#> [[4]]
#>            [,1]         [,2]        [,3]
#> [1,] -0.3819720 -0.001887537  0.00877901
#> [2,]  0.7961426 -0.476538944 -0.02036617
#> [3,] -0.4641852  0.649482522 -0.73753803
#> 
#> [[5]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3002220 -0.01198098 -0.06223524
#> [2,]  0.7810856 -0.54389046 -0.00553865
#> [3,] -0.4282115  0.75181463 -0.71607038
#> 
#> [[6]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3747898  0.02049325 -0.02505078
#> [2,]  0.8015257 -0.52168363 -0.00399670
#> [3,] -0.4478964  0.73241812 -0.66131745
#> 
#> [[7]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.4173450  0.03729497 -0.008754682
#> [2,]  0.8010541 -0.51885116  0.005999648
#> [3,] -0.4151946  0.70406628 -0.727095278
#> 
#> [[8]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3524343 -0.04499573 -0.021294682
#> [2,]  0.7735147 -0.51862260  0.002226638
#> [3,] -0.4122668  0.72039967 -0.695546670
#> 
#> [[9]]
#>            [,1]        [,2]         [,3]
#> [1,] -0.3568286  0.04395459 -0.009327648
#> [2,]  0.7308037 -0.49141152 -0.036315470
#> [3,] -0.4901935  0.74095200 -0.696247008
#> 
#> [[10]]
#>            [,1]        [,2]        [,3]
#> [1,] -0.3922283 -0.01480213 -0.01350045
#> [2,]  0.7497711 -0.51793784 -0.02178019
#> [3,] -0.4846939  0.75481409 -0.69685211
#> 
```
