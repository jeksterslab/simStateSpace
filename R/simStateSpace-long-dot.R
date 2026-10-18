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
  if (!eta) {
    varnames <- varnames[!(varnames %in% eta_names)]
  }
  # fill a preallocated matrix; binding per-person blocks copies the data
  # once per person
  rows <- vapply(
    X = x$data,
    FUN = function(i) {
      length(i$time)
    },
    FUN.VALUE = numeric(1)
  )
  last <- cumsum(rows)
  out <- matrix(
    data = 0.0,
    nrow = last[length(last)],
    ncol = length(varnames),
    dimnames = list(
      NULL,
      varnames
    )
  )
  cols_y <- 2 + seq_len(k)
  cols_eta <- 2 + k + seq_len(p)
  cols_x <- 2 + k + eta * p + seq_len(j)
  for (i in seq_along(x$data)) {
    data_i <- x$data[[i]]
    index <- (last[i] - rows[i] + 1):last[i]
    out[index, 1] <- data_i$id
    out[index, 2] <- data_i$time
    out[index, cols_y] <- data_i$y
    if (eta) {
      out[index, cols_eta] <- data_i$eta
    }
    if (j > 0) {
      out[index, cols_x] <- data_i$x
    }
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
