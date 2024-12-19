.DynrData <- function(object) {
  data <- .Long(
    x = object,
    eta = FALSE
  )
  k <- attributes(data)$k
  j <- attributes(data)$j
  if (j > 0) {
    output <- dynr::dynr.data(
      dataframe = as.data.frame(
        x = data
      ),
      id = "id",
      time = "time",
      observed = paste0("y", seq_len(k)),
      covariates = paste0("x", seq_len(j))
    )
  } else {
    output <- dynr::dynr.data(
      dataframe = as.data.frame(
        x = data
      ),
      id = "id",
      time = "time",
      observed = paste0("y", seq_len(k))
    )
  }
  return(output)
}
