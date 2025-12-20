# Spectral Abscissa

Returns the maximum real part of the eigenvalues of a square matrix. For
continuous-time stability (Hurwitz), a matrix is stable if the spectral
abscissa is strictly less than 0.

## Usage

``` r
SpectralAbscissa(x)
```

## Arguments

- x:

  Numeric square matrix.

## Value

Numeric value \\\alpha(x) = \max \Re(\lambda_i(x))\\.

## Author

Ivan Jacob Agaloos Pesigan

## Examples

``` r
# Hurwitz-stable (spectral abscissa < 0):
x <- matrix(
  data = c(
    -0.5, -0.2,
     1.0, -0.3
  ),
  nrow = 2
)
SpectralAbscissa(x = x) # < 0
#> [1] -0.4

# Unstable (spectral abscissa > 0):
x <- matrix(
  data = c(
     0.10, 0.50,
    -0.40, 0.20
  ),
  nrow = 2
)
SpectralAbscissa(x = x) # > 0
#> [1] 0.15
```
