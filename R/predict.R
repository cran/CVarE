#' Predict method for CVE Fits.
#'
#' Predict response using projected data. The forward model \eqn{g(B' X)} is
#' estimated with \code{\link{mars}} in the \pkg{mda} package.
#'
#' @param object an object of class \code{"cve"}, usually, a result of a call to
#'          \code{\link{cve}} or \code{\link{cve.call}}.
#' @param newdata Matrix of new predictor values, \eqn{C}.
#' @param k dimension of SDR space to be used for data projection.
#' @param ... further arguments passed to \code{\link{mars}}.
#'
#' @return prediced respone(s) for \code{newdata}.
#'
#' @examples
#' # create B for simulation
#' B <- rep(1, 5) / sqrt(5)
#'
#' set.seed(21)
#' # creat predictor data x ~ N(0, I_p)
#' x <- matrix(rnorm(500), 100)
#'
#' # simulate response variable
#' #    y = f(B'x) + err
#' # with f(x1) = x1 and err ~ N(0, 0.25^2)
#' y <- x %*% B + 0.25 * rnorm(100)
#'
#' x.train <- x[1:80, ]
#' x.test  <- x[81:100, ]
#' y.train <- y[1:80, ]
#' y.test  <- y[81:100, ]
#'
#' # calculate cve with method 'simple' for k = 1
#' cve.obj.simple <- cve(y.train ~ x.train, k = 1) 
#'
#' # predict y.test from x.test
#' yhat <- predict(cve.obj.simple, x.test, 1)
#'
#' # plot prediction against y.test
#' plot(yhat, y.test)
#' @seealso \code{\link{cve}}, \code{\link{cve.call}} and \pkg{\link{mars}}.
#'
#' @rdname predict.cve
#'
#' @importFrom mda mars
#' @method predict cve
#' @export
predict.cve <- function(object, newdata, k, ...) {
    if (missing(newdata)) {
        stop("No data supplied.")
    }
    if (missing(k)) {
        stop("No dimension supplied.")
    }
    
    if (!is.matrix(newdata)) {
        newdata <- matrix(newdata, nrow = 1L)
    }

    B <- object$res[[as.character(k)]]$B

    model <- mda::mars(object$X %*% B, object$Y)
    predict(model, newdata %*% B)
}
