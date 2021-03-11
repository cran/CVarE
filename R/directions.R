#' @export
directions <- function(object, k, ...) {
    UseMethod("directions")
}

#' Computes projected training data \code{X} for given dimension `k`.
#'
#' Returns \eqn{B'X}. That is, it computes the projection of the \eqn{n x p}
#' design matrix \eqn{X} on the column space of \eqn{B} of dimension \eqn{k}.
#'
#' @param object an object of class \code{"cve"}, usually, a result of a call to
#'          \code{\link{cve}} or \code{\link{cve.call}}.
#' @param k SDR dimension to use for projection.
#' @param ... ignored (no additional arguments).
#'
#' @return the \eqn{n\times k}{n x k} dimensional matrix \eqn{X B} where \eqn{B}
#'  is the  cve-estimate for dimension \eqn{k}.
#'
#' @examples
#' # create B for simulation (k = 1)
#' B <- rep(1, 5) / sqrt(5)
#' set.seed(21)
#' # creat predictor data x ~ N(0, I_p)
#' x <- matrix(rnorm(500), 100, 5)
#' # simulate response variable
#' #    y = f(B'x) + err
#' # with f(x1) = x1 and err ~ N(0, 0.25^2)
#' y <- x %*% B + 0.25 * rnorm(100)
#' # calculate cve with method 'mean' for k = 1
#' set.seed(21)
#' cve.obj.mean <- cve(y ~ x, k = 1, method = 'mean')
#' # get projected data for k = 1
#' x.proj <- directions(cve.obj.mean, k = 1)
#' # plot y against projected data
#' plot(x.proj, y)
#'
#' @seealso \code{\link{cve}}
#'
#' @method directions cve
#' @aliases directions directions.cve
#' @export
directions.cve <- function(object, k, ...) {
    if (!(k %in% names(object$res))) {
        stop("SDR directions for requested dimension `k` not computed.")
    }
    return(object$X %*% object$res[[as.character(k)]]$B)
}
