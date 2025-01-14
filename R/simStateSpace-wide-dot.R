.Wide <- function(x,
                  eta) {
  long <- .Long(
    x = x,
    eta = eta
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
  attributes(out)$k <- dims$k
  attributes(out)$p <- dims$p
  attributes(out)$j <- dims$j
  return(out)
}
