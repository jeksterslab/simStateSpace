.FormulaLinSDE <- function(p,
                           intercept) {
  formula <- lapply(
    X = seq_len(p),
    FUN = function(i) {
      terms <- paste0(
        "(",
        "phi_",
        i,
        "_",
        seq_len(p),
        " * eta_",
        seq_len(p),
        ")",
        collapse = " + "
      )
      return(
        paste0(
          "eta_",
          i,
          " ~ ",
          terms
        )
      )
    }
  )
  if (intercept) {
    return(
      lapply(
        X = seq_len(length(formula)),
        FUN = function(i) {
          paste0(
            formula[[i]],
            " + ",
            "iota_",
            i,
            "_1"
          )
        }
      )
    )
  } else {
    return(
      formula
    )
  }
}
