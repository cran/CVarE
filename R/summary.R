#' Prints summary statistics of the \eqn{L} \code{cve} component.
#'
#' Prints a summary statistics of the \code{L} component of a \code{cve} object #' for \code{k = min.dim, ..., max.dim}.
#'
#' @param object an object of class \code{"cve"}, usually, a result of a call to
#'          \code{\link{cve}} or \code{\link{cve.call}}.
#' @param ... ignored.
#'
#' @return No return value, prints human readable summary.
#'
#' @examples
#' # create B for simulation
#' B <- rep(1, 5) / sqrt(5)
#'
#' set.seed(21)
#' # create predictor data x ~ N(0, I_p)
#' x <- matrix(rnorm(500), 100)
#'
#' # simulate response variable
#' #    y = f(B'x) + err
#' # with f(x1) = x1 and err ~ N(0, 0.25^2)
#' y <- x %*% B + 0.25 * rnorm(100)
#' 
#' # calculate cve for unknown reduction dimension.
#' cve.obj.simple <- cve(y ~ x)
#' 
#' summary(cve.obj.simple)
#'
#' @method summary cve
#' @export
summary.cve <- function(object, ...) {
    cat('Summary of CVE result - Method: "', object$method, '"\n',
        '\n',
        'Dataset size:   ', nrow(object$X), '\n',
        'Data Dimension: ', ncol(object$X), '\n',
        # 'SDR Dimension:  ', object$k, '\n',
        # 'loss:           ', object$loss, '\n',
        '\n',
        'Called via:\n',
        '    ',
        sep='')
    print(object$call)

    L <- c()
    k <- c()
    for (dr.k in object$res) {
        if (class(dr.k) == 'cve.k') {
            k <- c(k, as.character(dr.k$k))
            L <- c(L, dr.k$L)
        }
    }
    L <- matrix(L, ncol = length(k))
    S <- apply(L, 2, summary)
    colnames(S) <- k
    cat('\n')
    print(S)
}
