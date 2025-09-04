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
  attributes(out)$n <- dims$n
  attributes(out)$m <- dims$m
  attributes(out)$k <- dims$k
  attributes(out)$p <- dims$p
  attributes(out)$j <- dims$j
  out
}
