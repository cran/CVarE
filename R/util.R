#' Random sample from Stiefel manifold.
#'
#' Draws a random sample from the invariant measure on the Stiefel manifold
#' \eqn{S(p, q)}.
#'
#' @param p row dimension
#' @param q col dimension
#' @return A \eqn{p \times q}{p x q} semi-orthogonal matrix.
#' @examples
#' V <- rStiefel(6, 4)
#' @export
rStiefel <- function(p, q) {
    return(qr.Q(qr(matrix(rnorm(p * q, 0, 1), p, q))))
}

#' Null space basis of given matrix `V`
#'
#' @param V `(p, q)` matrix
#' @return Semi-orthogonal `(p, p - q)` matrix spaning the null space of `V`.
#' @keywords internal
#' @export
null <- function(V) {
    tmp <- qr(V)
    set <- if(tmp$rank == 0L) seq_len(ncol(V)) else -seq_len(tmp$rank)
    return(qr.Q(tmp, complete = TRUE)[, set, drop = FALSE])
}
