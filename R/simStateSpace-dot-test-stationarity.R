.TestStationarity <- function(A) {
  eigen <- eigen(x = A)
  return(
    all(
      abs(eigen$values) < 1
    )
  )
}
