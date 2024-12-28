.Vec2Mat <- function(x) {
  prefix <- unique(
    gsub("_\\d+_\\d+$", "", names(x))
  )
  output <- lapply(
    X = prefix,
    FUN = function(prefix) {
      if (
        prefix %in% c(
          "theta",
          "psi",
          "sigma",
          "sigma0"
        )
      ) {
        return(
          .Vec2SymMat(
            x = x,
            prefix = prefix,
            row_names = NULL,
            col_names = NULL
          )
        )
      } else {
        return(
          .Vec2FullMat(
            x = x,
            prefix = prefix,
            row_names = NULL,
            col_names = NULL
          )
        )
      }
    }
  )
  names(output) <- prefix
  return(output)
}
