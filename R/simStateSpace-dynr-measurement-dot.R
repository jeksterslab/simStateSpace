.DynrMeasurement <- function(lambda,
                             nu) {
  p <- dim(lambda)[2]
  k <- dim(lambda)[1]
  values_int <- matrix(
    data = nu,
    ncol = 1
  )
  params_int <- matrix(
    data = paste0(
      "nu_",
      seq_len(k),
      "_1"
    ),
    ncol = 1
  )
  for (i in seq_len(dim(params_int)[1])) {
    if (values_int[i, 1] == 0) {
      params_int[i, 1] <- "fixed"
    }
  }
  return(
    dynr::prep.measurement(
      values.load = lambda,
      params.load = .LambdaLabel(lambda = lambda),
      state.names = paste0("eta_", seq_len(p)),
      obs.names = paste0("y", seq_len(k)),
      values.int = values_int,
      params.int = params_int
    )
  )
}
