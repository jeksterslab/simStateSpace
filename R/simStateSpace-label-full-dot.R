.LabelFull <- function(p,
                       label) {
  chars <- seq_len(p)
  return(
    outer(
      X = chars,
      Y = chars,
      FUN = function(x, y) {
        return(
          paste0(
            label,
            "_",
            x,
            "_",
            y
          )
        )
      }
    )
  )
}
