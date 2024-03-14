#' @inheritParams SimSSMFixed
Mu <- function(alpha,
               beta,
               nu) {
  .Mu0(
    alpha = alpha,
    beta = beta,
    nu = nu
  )
}
