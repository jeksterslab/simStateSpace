.FormulaOU <- function(p,
                       intercept) {
  if (intercept) {
    return(
      lapply(
        X = seq_len(p),
        FUN = function(i) {
          terms <- paste0(
            "(phi_",
            i,
            "_",
            seq_len(p),
            " * (eta_",
            seq_len(p),
            " - mu_",
            seq_len(p),
            "_1",
            "))",
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
    )
  } else {
    return(
      lapply(
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
    )
  }
}
