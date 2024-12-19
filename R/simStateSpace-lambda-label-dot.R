.LambdaLabel <- function(lambda) {
  p <- dim(lambda)[2]
  k <- dim(lambda)[1]
  lambda_labels <- matrix(
    data = "",
    nrow = k,
    ncol = p
  )
  for (j in seq_len(p)) {
    for (i in seq_len(k)) {
      lambda_labels[i, j] <- paste0(
        "lambda_",
        i,
        j
      )
    }
  }
  # set "fixed" for 0s in lambda
  for (j in seq_len(p)) {
    for (i in seq_len(k)) {
      if (lambda[i, j] == 0) {
        lambda_labels[i, j] <- "fixed"
      }
    }
  }
  # identify the anchor variable (the row with the value 1 that does not load on other columns)
  for (j in seq_len(p)) {
    candidate_rows <- which(lambda[, j] == 1) # Rows with 1 in column j
    for (i in candidate_rows) {
      # Check if the row does not load on other columns
      if (all(lambda[i, -j] == 0)) {
        lambda_labels[i, j] <- "fixed" # Mark this as the anchor
        break # Stop after finding the first valid anchor for this column
      }
    }
  }
  return(lambda_labels)
}
