# Changelog

## simStateSpace 1.2.15

### Patch

- Added the
  [`SSMInterceptEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMInterceptEta.md)
  and
  [`SSMInterceptY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMInterceptY.md)
  functions.

## simStateSpace 1.2.14

CRAN release: 2026-01-10

### Patch

- `CXX_STD = CXX17` in Makevars

## simStateSpace 1.2.13

### Patch

- Added the
  [`SimBetaNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaNCovariate.md)
  and
  [`SimPhiNCovariate()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiNCovariate.md)
  functions.

## simStateSpace 1.2.12

CRAN release: 2025-10-10

### Patch

- Added citation to Pesigan, I. J. A., Russell, M. A., & Chow, S.-M.
  (2025). Inferences and effect sizes for direct, indirect, and total
  effects in continuous-time mediation models. *Psychological Methods*.
  <https://doi.org/10.1037/met0000779>.

## simStateSpace 1.2.11

CRAN release: 2025-09-26

### Patch

- Added the
  [`SimAlphaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimAlphaN.md),
  [`SimIotaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimIotaN.md),
  and
  [`SimNuN()`](https://github.com/jeksterslab/simStateSpace/reference/SimNuN.md)
  functions.
- Added the
  [`LinSDECovEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDECovEta.md),
  [`LinSDECovY()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDECovY.md),
  [`LinSDEMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEMeanEta.md)
  and
  [`LinSDEMeanY()`](https://github.com/jeksterslab/simStateSpace/reference/LinSDEMeanY.md)
  functions. Removed the `LinSDECov()` and `LinSDEMean()` functions.
- Added the
  [`SSMCovEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMCovEta.md),
  [`SSMCovY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMCovY.md),
  [`SSMMeanEta()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanEta.md)
  and
  [`SSMMeanY()`](https://github.com/jeksterslab/simStateSpace/reference/SSMMeanY.md)
  functions.
- Added the
  [`SimCovN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovN.md)
  and
  [`SimCovDiagN()`](https://github.com/jeksterslab/simStateSpace/reference/SimCovDiagN.md)
  functions.
- Added the
  [`SimBetaN2()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN2.md)
  and
  [`SimPhiN2()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN2.md)
  functions.
- Transition matrices from
  [`SimBetaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN.md)
  can now be generated with optional bounds.
- Drift matrices from
  [`SimPhiN()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN.md)
  can now be generated with optional bounds.

## simStateSpace 1.2.10

CRAN release: 2025-03-30

### Patch

- Make sure the result of `LinSDECov()` is symmetric.

## simStateSpace 1.2.9

CRAN release: 2025-02-14

### Patch

- Added the `LinSDECov()` and `LinSDEMean()` functions.

## simStateSpace 1.2.8

CRAN release: 2025-01-23

### Patch

- Moved bootstrap components to a dedicated package, `bootStateSpace`,
  for better modularity.

## simStateSpace 1.2.7

CRAN release: 2025-01-08

### Patch

- Added the `PBSSMFixed()`, `PBSSMVARFixed()`, `PBSSMLinSDEFixed()`, and
  `PBSSMOUFixed()` functions.
- Added `SystemRequirements: GSL (>= 2.6)` in the DESCRIPTION file for
  the `dynr` package.

## simStateSpace 1.2.3

CRAN release: 2024-11-26

### Patch

- Replaced `arma::inv` with `arma::solve`.

## simStateSpace 1.2.2

CRAN release: 2024-06-21

### Patch

- Added the
  [`TestStationarity()`](https://github.com/jeksterslab/simStateSpace/reference/TestStationarity.md),
  [`TestStability()`](https://github.com/jeksterslab/simStateSpace/reference/TestStability.md),
  and
  [`TestPhi()`](https://github.com/jeksterslab/simStateSpace/reference/TestPhi.md)
  functions.
- Added the
  [`SimBetaN()`](https://github.com/jeksterslab/simStateSpace/reference/SimBetaN.md)
  and
  [`SimPhiN()`](https://github.com/jeksterslab/simStateSpace/reference/SimPhiN.md)
  functions.
- The
  [`SimSSMLinSDEIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEIVary.md)
  and
  [`SimSSMOUIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUIVary.md)
  functions can now accept a matrix of zeros for the argument `sigma_l`.
- The
  [`SimSSMLinSDEIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEIVary.md)
  function can now accept a vector of zeros for the argument `iota`.
- The
  [`SimSSMOUIVary()`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMOUIVary.md)
  function can now accept a vector of zeros for the argument `mu`.

## simStateSpace 1.2.1

CRAN release: 2024-05-14

### Patch

- The `LinSDE2SSM` function can now accept a matrix of zeros for the
  argument `sigma_l`.
- The `LinSDE2SSM` function can now accept a vector of zeros for the
  argument `iota`.

## simStateSpace 1.2.0

CRAN release: 2024-02-16

### Patch

- Added functions to generate linear stochastic differential equation
  models.

## simStateSpace 1.1.0

CRAN release: 2024-01-15

### Minor

- Added functions to generate data for models with covariates.
- Added functions to generate data for linear growth curve models.
- Added `print`, `as.data.frame`, `as.matrix`, and `plot` methods.

## simStateSpace 1.0.1

CRAN release: 2023-11-17

### Patch

- Updates to package documentation.

## simStateSpace 1.0.0

- Initial release.
