#' Extracts estimated SDR basis.
#'
#' Returns the SDR basis matrix for dimension \code{k}, i.e. returns the
#' cve-estimate of \eqn{B} with dimension \eqn{p\times k}{p x k}.
#'
#' @param object an object of class \code{"cve"}, usually, a result of a call to
#'          \code{\link{cve}} or \code{\link{cve.call}}.
#' @param k the SDR dimension.
#' @param ... ignored (no additional arguments).
#'
#' @return The matrix \eqn{B} of dimensions \eqn{p\times k}{p x k}.
#'
#' @examples
#' # set dimensions for simulation model
#' p <- 8 # sample dimension
#' k <- 2 # real dimension of SDR subspace
#' n <- 100 # samplesize
#' # create B for simulation
#' b1 <- rep(1 / sqrt(p), p)
#' b2 <- (-1)^seq(1, p) / sqrt(p)
#' B <- cbind(b1, b2)
#'
#' set.seed(21)
#' # creat predictor data x ~ N(0, I_p)
#' x <- matrix(rnorm(n * p), n, p)
#' # simulate response variable
#' #    y = f(B'x) + err
#' # with f(x1, x2) = x1^2 + 2 * x2 and err ~ N(0, 0.1^2)
#' y <- (x %*% b1)^2 + 2 * (x %*% b2) + 0.1 * rnorm(100)
#' # calculate cve for k = 2, 3
#' cve.obj <- cve(y ~ x, min.dim = 2, max.dim = 3)
#' # get cve-estimate for B with dimensions (p, k = 2)
#' B2 <- coef(cve.obj, k = 2)
#'
#' # Projection matrix on span(B)
#' # equivalent to `B %*% t(B)` since B is semi-orthonormal
#' PB <- B %*% solve(t(B) %*% B) %*% t(B)
#' # Projection matrix on span(B2)
#' # equivalent to `B2 %*% t(B2)` since B2 is semi-orthonormal
#' PB2 <- B2 %*% solve(t(B2) %*% B2) %*% t(B2)
#' # compare estimation accuracy by Frobenius norm of difference of projections
#' norm(PB - PB2, type = 'F')
#'
#' @method coef cve
#' @aliases coef.cve
#' @rdname coef.cve
#' @export
coef.cve <- function(object, k, ...) {
    if (missing(k)) {
        Bs <- list()
        for (k in names(object$res)) {
            Bs[[k]] <- object$res[[k]]$B
        }
        return(Bs)
    } else if (k %in% names(object$res)) {
        return(object$res[[as.character(k)]]$B)
    } else {
        stop("Requested dimension `k` not computed.")
    }
}
