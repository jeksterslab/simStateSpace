.Long <- function(x,
                  eta) {
  first <- x$data[[1]]
  obs <- first$y
  k <- dim(obs)[2]
  y_names <- paste0("y", seq_len(k))
  varnames <- c(
    "id",
    "time",
    y_names
  )
  lats <- first$eta
  p <- dim(lats)[2]
  eta_names <- paste0("eta", seq_len(p))
  varnames <- c(
    varnames,
    eta_names
  )
  if (x$model$covariates) {
    covs <- first$x
    j <- dim(covs)[2]
    x_names <- paste0("x", seq_len(j))
    varnames <- c(
      varnames,
      x_names
    )
  } else {
    j <- 0
  }
  out <- lapply(
    X = x$data,
    FUN = function(x) {
      return(
        do.call(
          what = "cbind",
          args = x
        )
      )
    }
  )
  out <- do.call(
    what = "rbind",
    args = out
  )
  colnames(out) <- varnames
  if (!eta) {
    varnames <- varnames[!(varnames %in% eta_names)]
    out <- out[, varnames, drop = FALSE]
  }
  attributes(out)$n <- length(
    unique(out[, "id"])
  )
  attributes(out)$k <- k
  attributes(out)$p <- p
  attributes(out)$j <- j
  return(out)
}
