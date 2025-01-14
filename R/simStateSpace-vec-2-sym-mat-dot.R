.Vec2SymMat <- function(x,
                        prefix,
                        row_names,
                        col_names) {
  y <- x[
    grep(
      paste0(
        "^",
        prefix,
        "_"
      ),
      names(x)
    )
  ]
  index <- do.call(
    rbind,
    strsplit(
      names(y),
      "_"
    )
  )
  row <- as.numeric(index[, 2])
  col <- as.numeric(index[, 3])
  y_matrix <- matrix(
    data = 0,
    nrow = max(row),
    ncol = max(col)
  )
  for (i in seq_len(length(row))) {
    y_matrix[row[i], col[i]] <- y[i]
    y_matrix[col[i], row[i]] <- y[i]
  }
  if (!is.null(row_names)) {
    rownames(y_matrix) <- row_names
  }
  if (!is.null(col_names)) {
    colnames(y_matrix) <- col_names
  }
  return(y_matrix)
}
