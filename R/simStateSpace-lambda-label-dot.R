.LambdaLabel <- function(lambda) {
  p <- dim(lambda)[2]
  k <- dim(lambda)[1]
  .GenerateLabels <- function(k,
                              p) {
    outer(
      X = seq_len(k),
      Y = seq_len(p),
      FUN = Vectorize(
        FUN = function(i, j) {
          paste0("lambda_", i, j)
        }
      )
    )
  }
  .SetFixed <- function(labels,
                        lambda) {
    zero_indices <- which(
      x = lambda == 0,
      arr.ind = TRUE
    )
    labels[zero_indices] <- "fixed"
    return(labels)
  }
  .SetAnchors <- function(labels,
                          lambda,
                          p) {
    for (j in seq_len(p)) {
      candidate_rows <- which(lambda[, j] == 1)
      for (i in candidate_rows) {
        if (all(lambda[i, -j] == 0)) {
          labels[i, j] <- "fixed"
          break
        }
      }
    }
    return(labels)
  }
  labels <- .GenerateLabels(
    k = k,
    p = p
  )
  labels <- .SetFixed(
    labels = labels,
    lambda = lambda
  )
  return(
    .SetAnchors(
      labels = labels,
      lambda = lambda,
      p = p
    )
  )
}
