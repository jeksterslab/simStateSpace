#' Print Method for an Object of Class
#' `statespacepb`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Prints a matrix of
#'   estimates,
#'   standard errors,
#'   number of bootstrap replications,
#'   and
#'   confidence intervals.
#'
#' @param x Object of Class `statespacepb`.
#' @param alpha Numeric vector.
#'   Significance level \eqn{\alpha}.
#'   If `alpha = NULL`,
#'   use the argument `alpha` used in `x`.
#' @inheritParams summary.statespacepb
#'
#' @keywords methods
#' @export
print.statespacepb <- function(x,
                               alpha = NULL,
                               type = "pc",
                               digits = 4,
                               ...) {
  cat("Call:\n")
  base::print(x$call)
  cat(
    paste0(
      "\n",
      "Parametric bootstrap confidence intervals.",
      "\n",
      "type = ",
      "\"",
      type,
      "\"",
      "\n"
    )
  )
  base::print(
    round(
      .PBCI(
        object = x,
        alpha = alpha,
        type = type
      ),
      digits = digits
    )
  )
}

#' Summary Method for an Object of Class
#' `statespacepb`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns a matrix of
#'   estimates,
#'   standard errors,
#'   number of bootstrap replications,
#'   and
#'   confidence intervals.
#'
#' @param object Object of Class `statespacepb`.
#' @param ... additional arguments.
#' @param alpha Numeric vector.
#'   Significance level \eqn{\alpha}.
#'   If `alpha = NULL`,
#'   use the argument `alpha` used in `object`.
#' @param type Charater string.
#'   Confidence interval type, that is,
#'   `type = "pc"` for percentile;
#'   `type = "bc"` for bias corrected.
#' @param digits Digits to print.
#'
#' @keywords methods
#' @export
summary.statespacepb <- function(object,
                                 alpha = NULL,
                                 type = "pc",
                                 digits = 4,
                                 ...) {
  cat("Call:\n")
  base::print(object$call)
  cat(
    paste0(
      "\n",
      "Parametric bootstrap confidence intervals.",
      "\n",
      "type = ",
      "\"",
      type,
      "\"",
      "\n"
    )
  )
  return(
    round(
      .PBCI(
        object = object,
        alpha = alpha,
        type = type
      ),
      digits = digits
    )
  )
}

#' Sampling Variance-Covariance Matrix Method for an Object of Class
#' `statespacepb`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns the variance-covariance matrix of estimates.
#'
#' @inheritParams summary.statespacepb
#'
#' @keywords methods
#' @export
vcov.statespacepb <- function(object,
                              ...) {
  return(
    object$vcov
  )
}

#' Estimated Parameter Method for an Object of Class
#' `statespacepb`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns a vector of estimated parameters.
#'
#' @inheritParams summary.statespacepb
#'
#' @keywords methods
#' @export
coef.statespacepb <- function(object,
                              ...) {
  return(
    object$est
  )
}

#' Confidence Intervals Method for an Object of Class
#' `statespacepb`
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns a matrix of confidence intervals.
#'
#' @inheritParams summary.statespacepb
#' @param parm a specification of which parameters
#'   are to be given confidence intervals,
#'   either a vector of numbers or a vector of names.
#'   If missing, all parameters are considered.
#' @param level the confidence level required.
#'
#' @keywords methods
#' @export
confint.statespacepb <- function(object,
                                 parm = NULL,
                                 level = 0.95,
                                 type = "pc",
                                 ...) {
  if (is.null(parm)) {
    parm <- seq_len(
      length(
        object$est
      )
    )
  }
  ci <- .PBCI(
    object = object,
    alpha = 1 - level[1],
    type = type
  )[parm, 4:5, drop = FALSE]
  varnames <- colnames(ci)
  varnames <- gsub(
    pattern = "%",
    replacement = " %",
    x = varnames
  )
  colnames(ci) <- varnames
  return(
    ci
  )
}
