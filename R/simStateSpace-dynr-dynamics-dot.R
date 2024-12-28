.DynrDynamics <- function(dynamics,
                          intercept,
                          continuous,
                          ou) {
  if (all(intercept == 0)) {
    intercept <- FALSE
  } else {
    intercept_value <- intercept
    intercept <- TRUE
  }
  if (continuous) {
    # phi
    dynamics_label <- .LabelFull(
      p = dim(dynamics)[1],
      label = "phi"
    )
    if (ou) {
      formula <- .FormulaOU(
        p = dim(dynamics)[1],
        intercept = intercept
      )
    } else {
      formula <- .FormulaLinSDE(
        p = dim(dynamics)[1],
        intercept = intercept
      )
    }
  } else {
    # beta
    dynamics_label <- .LabelFull(
      p = dim(dynamics)[1],
      label = "beta"
    )
    formula <- .FormulaVAR(
      p = dim(dynamics)[1],
      intercept = intercept
    )
  }
  startval <- .Vec(dynamics)
  names(startval) <- .Vec(dynamics_label)
  if (intercept) {
    if (continuous) {
      if (ou) {
        intercept_label <- paste0(
          "mu_",
          seq_len(
            length(intercept_value)
          ),
          "_1"
        )
      } else {
        intercept_label <- paste0(
          "iota_",
          seq_len(
            length(intercept_value)
          ),
          "_1"
        )
      }
    } else {
      intercept_label <- paste0(
        "alpha_",
        seq_len(
          length(intercept_value)
        ),
        "_1"
      )
    }
    names(intercept_value) <- intercept_label
    startval <- c(
      startval,
      intercept_value
    )
  }
  return(
    list(
      dynamics = dynr::prep.formulaDynamics(
        formula = lapply(
          X = formula,
          FUN = stats::as.formula
        ),
        startval = startval,
        isContinuousTime = continuous
      ),
      startval = startval
    )
  )
}
