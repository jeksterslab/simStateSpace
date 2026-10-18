.CheckCTVAR <- function(iota, phi, sigma_l, delta_t) {
  tol <- sqrt(.Machine$double.eps)
  if (is.matrix(iota)) {
    dim(iota) <- NULL
  }
  stopifnot(
    "`iota` must be a vector." = is.null(dim(iota)),
    "`iota` must be numeric." = is.numeric(iota),
    "`iota` must contain only finite values." = all(is.finite(iota)),
    "`phi` must be a list." = is.list(phi),
    "`phi` must have positive length." = length(phi) > 0L,
    "`phi` must be a list of matrices." = all(
      vapply(
        X = phi,
        FUN = is.matrix,
        FUN.VALUE = logical(1)
      )
    ),
    "`phi` must be a list of numeric matrices." = all(
      vapply(
        X = phi,
        FUN = is.numeric,
        FUN.VALUE = logical(1)
      )
    ),
    "`phi` must contain only finite values." = all(
      vapply(
        X = phi,
        FUN = function(x) {
          all(is.finite(x))
        },
        FUN.VALUE = logical(1)
      )
    ),
    "`sigma_l` must be a matrix." = is.matrix(sigma_l),
    "`sigma_l` must be numeric." = is.numeric(sigma_l),
    "`sigma_l` must contain only finite values." = all(is.finite(sigma_l)),
    "`delta_t` must be numeric." = is.numeric(delta_t),
    "`delta_t` must have length 1." = length(delta_t) == 1L,
    "`delta_t` must be finite." = is.finite(delta_t),
    "`delta_t` must be positive." = delta_t > 0.0
  )
  q <- nrow(phi[[1]])
  stopifnot(
    "`phi[[1]]` must have positive dimensions." = q > 0L,
    "`phi[[1]]` must be square." = nrow(phi[[1]]) == ncol(phi[[1]]),
    "`iota` must have length equal to the number of CTVAR variables." =
      length(iota) == q,
    "All matrices in `phi` must have the same square dimensions." = all(
      vapply(
        X = phi,
        FUN = function(x) {
          identical(
            x = dim(x),
            y = c(q, q)
          )
        },
        FUN.VALUE = logical(1)
      )
    )
  )
  stopifnot(
    "`sigma_l` must have dimensions compatible with `iota` and `phi`." =
      identical(
        x = dim(sigma_l),
        y = c(q, q)
      ),
    "`sigma_l` must be lower triangular." = all(
      abs(sigma_l[upper.tri(sigma_l)]) <= tol
    ),
    "`sigma_l` must have non-negative diagonal values." = all(
      diag(sigma_l) >= -tol
    )
  )
  invisible(NULL)
}
