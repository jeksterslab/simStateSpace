.TestPositiveDefinite <- function(A,
                                  tol = 1e-06) {
  eigen <- eigen(x = A, symmetric = TRUE)
  return(
    all(
      eigen$values >= -tol * abs(
        eigen$values[1L]
      )
    )
  )
}
