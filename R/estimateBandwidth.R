#' Bandwidth estimation for CVE.
#'
#' If no bandwidth or function for calculating it is supplied, the CVE method
#' defaults to using the following formula (version 1)
#' \deqn{%
#'    h = \frac{2 tr(\Sigma)}{p} (1.2 n^{\frac{-1}{4 + k}})^2}{%
#'    h = (2 * tr(\Sigma) / p) * (1.2 * n^(-1 / (4 + k)))^2}
#' Alternative version 2 is used for dimension prediction which is given by
#'    \deqn{%
#'    h = \frac{2 tr(\Sigma)}{p} \chi_k^{-1}(\frac{nObs - 1}{n - 1})}{%
#'    h = (2 * tr(\Sigma) / p) * \chi_k^-1((nObs - 1) / (n - 1))}
#' with \eqn{n} the sample size, \eqn{p} the dimension of \eqn{X} and 
#' \eqn{\Sigma} is \eqn{(n - 1) / n} times the sample covariance matrix of
#' \eqn{X}.
#'
#' @param X the \eqn{n\times p}{n x p} matrix of predictor values.
#' @param k the SDR dimension.
#' @param nObs number of points in a slice, only for version 2.
#' @param version either \code{1} or \code{2}.
#'
#' @return Estimated bandwidth \code{h}.
#'
#' @examples
#' # set dimensions for simulation model
#' p <- 5; k <- 1
#' # create B for simulation
#' B <- rep(1, p) / sqrt(p)
#' # samplsize
#' n <- 100
#' set.seed(21)
#' #creat predictor data x ~ N(0, I_p)
#' x <- matrix(rnorm(n * p), n, p)
#' # simulate response variable
#' #     y = f(B'x) + err
#' # with f(x1) = x1 and err ~ N(0, 0.25^2)
#' y <- x %*% B + 0.25 * rnorm(100)
#' # calculate cve with method 'simple' for k = 1
#' set.seed(21)
#' cve.obj.simple <- cve(y ~ x, k = k)
#' print(estimate.bandwidth(x, k = k))
#' @export
estimate.bandwidth <- function (X, k, nObs, version = 1L) {
    n <- nrow(X)
    p <- ncol(X)
    if (version == 1) {
        X_centered <- scale(X, center = TRUE, scale = FALSE)
        Sigma <- crossprod(X_centered, X_centered)/n
        return((2 * sum(diag(Sigma))/p) * (1.2 * n^(-1/(4 + k)))^2)
    } else if (version == 2) {
        X_c <- scale(X, center = TRUE, scale = FALSE)
        return(2 * qchisq((nObs - 1) / (n - 1), k) * sum(X_c^2) / (n * p))
    } else {
        stop("Unknown version.")
    }
}
