#' Confidence Intervals
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @return Returns a matrix of
#'   estimates,
#'   standard errors,
#'   and
#'   confidence intervals.
#'
#' @param object Object of class `statespacepb`.
#' @param alpha Numeric vector.
#'   Significance level \eqn{\alpha}.
#' @param type Charater string.
#'   Confidence interval type, that is,
#'   `type = "pc"` for percentile;
#'   `type = "bc"` for bias corrected.
#'
#' @family Beta Nonparametric Bootstrap Functions
#' @keywords betaNB ci internal
#' @noRd
.PBCI <- function(object,
                  alpha = NULL,
                  type = "pc") {
  if (is.null(alpha)) {
    alpha <- object$args$alpha_level
  }
  probs <- .PCProbs(alpha = alpha)
  thetahatstar <- do.call(
    what = "rbind",
    args = object$thetahatstar
  )
  thetahatstar <- unname(thetahatstar)
  colnames(thetahatstar) <- names(object$est)
  thetahat <- object$est
  k <- dim(thetahatstar)[2]
  ci <- vector(
    mode = "list",
    length = k
  )
  if (type == "pc") {
    for (i in seq_len(k)) {
      ci[[i]] <- .PCCI(
        thetahatstar = thetahatstar[, i],
        thetahat = thetahat[[i]],
        probs = probs
      )
    }
  }
  if (type == "bc") {
    z1 <- .Z1(probs = probs)
    for (i in seq_len(k)) {
      ci[[i]] <- .BCCI(
        thetahatstar = thetahatstar[, i],
        thetahat = thetahat[[i]],
        probs = probs,
        z0 = .Z0(
          thetahatstar = thetahatstar[, i],
          thetahat = thetahat[[i]]
        ),
        z1 = z1
      )
    }
  }
  ci <- do.call(
    what = "rbind",
    args = ci
  )
  rownames(ci) <- names(thetahat)
  return(
    ci
  )
}
