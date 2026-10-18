.Wide <- function(x,
                  eta,
                  burnin = 0,
                  reset_time = TRUE) {
  long <- .Long(
    x = x,
    eta = eta,
    burnin = burnin,
    reset_time = reset_time
  )
  dims <- attributes(long)
  times <- unique(long[, "time"])
  m <- length(times)
  rows <- nrow(long) / m
  varying <- colnames(long)[-c(1, 2)]
  # every person has the same occasions in the same order, so each person's
  # block can be laid out as one row (same columns as stats::reshape())
  balanced <- rows == dims$n && rows == floor(rows) &&
    all(long[, "time"] == rep(x = times, times = rows)) &&
    all(long[, "id"] == rep(x = unique(long[, "id"]), each = m))
  if (balanced) {
    values <- long[, varying, drop = FALSE]
    dim(values) <- c(m, rows, length(varying))
    values <- aperm(
      a = values,
      perm = c(3, 1, 2)
    )
    dim(values) <- c(m * length(varying), rows)
    out <- cbind(
      unique(long[, "id"]),
      t(values)
    )
    dimnames(out) <- list(
      NULL,
      c(
        "id",
        paste(
          rep(x = varying, times = m),
          rep(x = times, each = length(varying)),
          sep = "_"
        )
      )
    )
  } else {
    out <- as.matrix(
      stats::reshape(
        data = as.data.frame(
          long
        ),
        timevar = "time",
        idvar = "id",
        direction = "wide",
        sep = "_"
      )
    )
    rownames(out) <- NULL
  }
  attributes(out)$n <- dims$n
  attributes(out)$m <- dims$m
  attributes(out)$k <- dims$k
  attributes(out)$p <- dims$p
  attributes(out)$j <- dims$j
  out
}
