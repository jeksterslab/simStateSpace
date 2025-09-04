.Long <- function(x,
                  eta,
                  burnin = 0,
                  reset_time = TRUE) {
  stopifnot(
    burnin >= 0
  )
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
      do.call(
        what = "cbind",
        args = x
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
  if (burnin > 0) {
    time <- sort(
      unique(
        out[, "time"]
      )
    )
    m <- length(time)
    if (burnin >= m) {
      stop(
        "`burnin` should not be greater than the measurement occasions.\n"
      )
    }
    out <- out[
      which(
        out[, "time"] %in% time[-seq_len(burnin)]
      ),
    ]
    if (reset_time) {
      out[, "time"] <- out[, "time"] - sort(
        unique(
          out[, "time"]
        )
      )[1]
    }
  }
  attributes(out)$n <- length(
    unique(out[, "id"])
  )
  attributes(out)$m <- length(
    unique(out[, "time"])
  )
  attributes(out)$k <- k
  attributes(out)$p <- p
  attributes(out)$j <- j
  out
}
