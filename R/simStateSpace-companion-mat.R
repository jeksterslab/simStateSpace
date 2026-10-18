#' Create a Block Companion Matrix
#'
#' Creates a block companion matrix from a list of coefficient matrices.
#'
#' @details
#' Let `x` be a list of coefficient matrices
#' \eqn{\{\Beta_1, \Beta_2, \ldots, \Beta_p\}},
#' where each \eqn{\Beta_j} is a \eqn{q \times q} matrix.
#' This function returns the \eqn{pq \times pq} block companion matrix
#' whose first block row contains the coefficient matrices and whose
#' lower block subdiagonal contains identity matrices.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param x A list of numeric coefficient matrices.
#'
#' @return A block companion matrix.
#'
#' @examples
#' beta <- list(
#'   matrix(
#'     data = c(
#'       0.40, 0.10,
#'       0.20, 0.30
#'     ),
#'     nrow = 2,
#'     ncol = 2
#'   ),
#'   matrix(
#'     data = c(
#'       0.10, 0.00,
#'       0.05, 0.10
#'     ),
#'     nrow = 2,
#'     ncol = 2
#'   )
#' )
#'
#' CompanionMat(x = beta)
#'
#' @family Simulation of State Space Models Data Functions
#' @keywords simStateSpace var
#' @export
CompanionMat <- function(x) {
  stopifnot(
    "`x` must be a list." = is.list(x),
    "`x` must not be a data frame." = !inherits(
      x = x,
      what = "data.frame"
    ),
    "`x` must have at least one element." = length(x) > 0L
  )

  is_matrix <- unlist(
    x = lapply(
      X = x,
      FUN = is.matrix
    ),
    use.names = FALSE
  )

  stopifnot(
    "Each element of `x` must be a matrix." = all(is_matrix)
  )

  is_numeric <- unlist(
    x = lapply(
      X = x,
      FUN = is.numeric
    ),
    use.names = FALSE
  )

  stopifnot(
    "Each element of `x` must be numeric." = all(is_numeric)
  )

  is_finite <- unlist(
    x = lapply(
      X = x,
      FUN = function(x_i) {
        all(is.finite(x_i))
      }
    ),
    use.names = FALSE
  )

  stopifnot(
    "Each element of `x` must contain only finite values." = all(is_finite)
  )

  is_square <- unlist(
    x = lapply(
      X = x,
      FUN = function(x_i) {
        nrow(x_i) == ncol(x_i)
      }
    ),
    use.names = FALSE
  )

  stopifnot(
    "Each element of `x` must be a square matrix." = all(is_square)
  )

  dim_x <- dim(x[[1]])

  stopifnot(
    "The matrices in `x` must have positive dimensions." = dim_x[1] > 0L
  )

  has_same_dim <- unlist(
    x = lapply(
      X = x,
      FUN = function(x_i) {
        identical(
          x = dim(x_i),
          y = dim_x
        )
      }
    ),
    use.names = FALSE
  )

  stopifnot(
    "All matrices in `x` must have the same dimensions." = all(has_same_dim)
  )

  x <- lapply(
    X = x,
    FUN = function(x_i) {
      storage.mode(x_i) <- "double"
      x_i
    }
  )

  .CompanionMat(x = x)
}
