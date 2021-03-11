#' Conditional Variance Estimator (CVE) Package.
#'
#' Conditional Variance Estimation (CVE) is a novel sufficient dimension
#' reduction (SDR) method for regressions satisfying \eqn{E(Y|X) = E(Y|B'X)},
#' where \eqn{B'X} is a lower dimensional projection of the predictors and
#' \eqn{Y} is a univariate response. CVE,
#' similarly to its main competitor, the mean average variance estimation
#' (MAVE), is not based on inverse regression, and does not require the
#' restrictive linearity and constant variance conditions of moment based SDR
#' methods. CVE is data-driven and applies to additive error regressions with
#' continuous predictors and link function. Let \eqn{X} be a real
#' \eqn{p}-dimensional covariate vector. We assume that the dependence of
#' \eqn{Y} and \eqn{X} is modelled by
#'
#' \deqn{Y = g(B'X) + \epsilon}
#'
#' where \eqn{X} is independent of \eqn{\epsilon} with positive definite
#' variance-covariance matrix \eqn{Var(X) = \Sigma_X}. \eqn{\epsilon} is a mean
#' zero random variable with finite \eqn{Var(\epsilon) = E(\epsilon^2)}, \eqn{g}
#' is an unknown, continuous non-constant function,
#' and \eqn{B = (b_1, ..., b_k)} is
#' a real \eqn{p \times k}{p x k} matrix of rank \eqn{k \leq p}{k <= p}. 
#' Without loss of generality \eqn{B} is assumed to be orthonormal.
#'
#' Further, the extended Ensemble Conditional Variance Estimation (ECVE) is
#' implemented which is a SDR method in regressions with continuous response and
#' predictors. ECVE applies to general non-additive error regression models.
#'
#' \deqn{Y = g(B'X, \epsilon)}
#'
#' It operates under the assumption that the predictors can be replaced by a
#' lower dimensional projection without loss of information.It is a
#' semiparametric forward regression model based exhaustive sufficient dimension
#' reduction estimation method that is shown to be consistent under mild
#' assumptions.
#'
#' @author Daniel Kapla, Lukas Fertl, Bura Efstathia
#'
#' @references
#'    [1] Fertl, L. and Bura, E. (2021) "Conditional Variance
#'          Estimation for Sufficient Dimension Reduction"
#'          <arXiv:2102.08782>
#'
#'    [2] Fertl, L. and Bura, E. (2021) "Ensemble Conditional Variance
#'          Estimation for Sufficient Dimension Reduction"
#'          <arXiv:2102.13435>
#'
#' @docType package
#' @useDynLib CVarE, .registration = TRUE
"_PACKAGE"

#' Conditional Variance Estimator (CVE).
#'
#' @description
#' This is the main function in the \code{CVE} package. It creates objects of
#' class \code{"cve"} to estimate the mean subspace. Helper functions that
#' require a \code{"cve"} object can then be applied to the output from this
#' function.
#'
#' Conditional Variance Estimation (CVE) is a sufficient dimension reduction
#' (SDR) method for regressions studying \eqn{E(Y|X)}, the conditional
#' expectation of a response \eqn{Y} given a set of predictors \eqn{X}. This
#' function provides methods for estimating the dimension and the subspace
#' spanned by the columns of a \eqn{p\times k}{p x k} matrix \eqn{B} of minimal
#' rank \eqn{k} such that
#'
#' \deqn{E(Y|X) = E(Y|B'X)}
#'
#' or, equivalently,
#' 
#' \deqn{Y = g(B'X) + \epsilon}
#'
#' where \eqn{X} is independent of \eqn{\epsilon} with positive definite
#' variance-covariance matrix \eqn{Var(X) = \Sigma_X}. \eqn{\epsilon} is a mean
#' zero random variable with finite \eqn{Var(\epsilon) = E(\epsilon^2)}, \eqn{g}
#' is an unknown, continuous non-constant function, and \eqn{B = (b_1,..., b_k)}
#' is a real \eqn{p \times k}{p x k} matrix of rank \eqn{k \leq p}{k <= p}.
#'
#' Both the dimension \eqn{k} and the subspace \eqn{span(B)} are unknown. The 
#' CVE method makes very few assumptions.
#'
#' A kernel matrix \eqn{\hat{B}}{Bhat} is estimated such that the column space
#' of \eqn{\hat{B}}{Bhat} should be close to the mean subspace \eqn{span(B)}.
#' The primary output from this method is a set of orthonormal vectors,
#' \eqn{\hat{B}}{Bhat}, whose span estimates \eqn{span(B)}.
#'
#' The method central implements the Ensemble Conditional Variance Estimation
#' (ECVE) as described in [2]. It augments the CVE method by applying an
#' ensemble of functions (parameter \code{func_list}) to the response to
#' estimate the central subspace. This corresponds to the generalization
#'
#' \deqn{F(Y|X) = F(Y|B'X)}
#'
#' or, equivalently,
#'
#' \deqn{Y = g(B'X, \epsilon)}
#'
#' where \eqn{F} is the conditional cumulative distribution function.
#'
#' @param formula an object of class \code{"formula"} which is a symbolic
#' description of the model to be fitted like \eqn{Y\sim X}{Y ~ X} where
#' \eqn{Y} is a \eqn{n}-dimensional vector of the response variable and
#' \eqn{X} is a \eqn{n\times p}{n x p} matrix of the predictors.
#' @param data an optional data frame, containing the data for the formula if
#' supplied like \code{data <- data.frame(Y, X)} with dimension
#' \eqn{n \times (p + 1)}{n x (p + 1)}. By default the variables are taken from
#' the environment from which \code{cve} is called.
#' @param method This character string specifies the method of fitting. The
#' options are
#' \itemize{
#'    \item \code{"mean"} method to estimate the mean subspace, see [1].
#'    \item \code{"central"} ensemble method to estimate the central subspace,
#'      see [2].
#'    \item \code{"weighted.mean"} variation of \code{"mean"} method with
#'      adaptive weighting of slices, see [1].
#'    \item \code{"weighted.central"} variation of \code{"central"} method with
#'      adaptive weighting of slices, see [2].
#' }
#' @param max.dim upper bounds for \code{k}, (ignored if \code{k} is supplied).
#' @param ... optional parameters passed on to \code{\link{cve.call}}.
#'
#' @return an S3 object of class \code{cve} with components:
#' \describe{
#'    \item{X}{design matrix of predictor vector used for calculating
#'      cve-estimate,}
#'    \item{Y}{\eqn{n}-dimensional vector of responses used for calculating
#'        cve-estimate,}
#'    \item{method}{Name of used method,}
#'    \item{call}{the matched call,}
#'    \item{res}{list of components \code{V, L, B, loss, h} for
#'       each \code{k = min.dim, ..., max.dim}. If \code{k} was supplied in the
#'       call \code{min.dim = max.dim = k}.
#'       \itemize{
#'           \item \code{B} is the cve-estimate with dimension
#'               \eqn{p\times k}{p x k}.
#'           \item \code{V} is the orthogonal complement of \eqn{B}.
#'           \item \code{L} is the loss for each sample seperatels such that
#'               it's mean is \code{loss}.
#'           \item \code{loss} is the value of the target function that is 
#'               minimized, evaluated at \eqn{V}.
#'           \item \code{h} bandwidth parameter used to calculate
#'               \code{B, V, loss, L}.
#'       }
#'    }
#' }
#'
#' @examples
#' # set dimensions for simulation model
#' p <- 5
#' k <- 2
#' # create B for simulation
#' b1 <- rep(1 / sqrt(p), p)
#' b2 <- (-1)^seq(1, p) / sqrt(p)
#' B <- cbind(b1, b2)
#' # sample size
#' n <- 100
#' set.seed(21)
#'
#' # creat predictor data x ~ N(0, I_p)
#' x <- matrix(rnorm(n * p), n, p)
#' # simulate response variable
#' #     y = f(B'x) + err
#' # with f(x1, x2) = x1^2 + 2 * x2 and err ~ N(0, 0.25^2)
#' y <- (x %*% b1)^2 + 2 * (x %*% b2) + 0.25 * rnorm(n)
#'
#' # calculate cve with method 'mean' for k unknown in 1, ..., 3
#' cve.obj.s <- cve(y ~ x, max.dim = 2) # default method 'mean'
#' # calculate cve with method 'weighed' for k = 2
#' cve.obj.w <- cve(y ~ x, k = 2, method = 'weighted.mean')
#' B2 <- coef(cve.obj.s, k = 2)
#'
#' # get projected X data (same as cve.obj.s$X %*% B2)
#' proj.X <- directions(cve.obj.s, k = 2)
#' #  plot y against projected data
#' plot(proj.X[, 1], y)
#' plot(proj.X[, 2], y)
#'
#' # creat 10 new x points and y according to model
#' x.new <- matrix(rnorm(10 * p), 10, p)
#' y.new <- (x.new %*% b1)^2 + 2 * (x.new %*% b2) + 0.25 * rnorm(10)
#' # predict y.new
#' yhat <- predict(cve.obj.s, x.new, 2)
#' plot(y.new, yhat)
#'
#' # projection matrix on span(B)
#' # same as B %*% t(B) since B is semi-orthogonal
#' PB <- B %*% solve(t(B) %*% B) %*% t(B)
#' # cve estimates for B with mean and weighted method
#' B.s <- coef(cve.obj.s, k = 2)
#' B.w <- coef(cve.obj.w, k = 2)
#' # same as B.s %*% t(B.s) since B.s is semi-orthogonal (same vor B.w)
#' PB.s <- B.s %*% solve(t(B.s) %*% B.s) %*% t(B.s)
#' PB.w <- B.w %*% solve(t(B.w) %*% B.w) %*% t(B.w)
#' # compare estimation accuracy of mean and weighted cve estimate by
#' # Frobenius norm of difference of projections.
#' norm(PB - PB.s, type = 'F')
#' norm(PB - PB.w, type = 'F')
#'
#' @seealso For a detailed description of \code{formula} see
#'      \code{\link{formula}}.
#'
#' @references
#'    [1] Fertl, L. and Bura, E. (2021) "Conditional Variance
#'          Estimation for Sufficient Dimension Reduction"
#'          <arXiv:2102.08782>
#'
#'    [2] Fertl, L. and Bura, E. (2021) "Ensemble Conditional Variance
#'          Estimation for Sufficient Dimension Reduction"
#'          <arXiv:2102.13435>
#'
#' @importFrom stats model.frame
#' @export
cve <- function(formula, data, method = "mean", max.dim = 10L, ...) {
    # check for type of `data` if supplied and set default
    if (missing(data)) {
        data <- environment(formula)
    } else if (!is.data.frame(data)) {
        stop("Parameter 'data' must be a 'data.frame' or missing.")
    }

    # extract `X`, `Y` from `formula` with `data`
    model <- stats::model.frame(formula, data)
    Y <- stats::model.response(model, "double")
    X <- stats::model.matrix(model, data)
    if ("(Intercept)" %in% colnames(X)) {
        X <- X[, "(Intercept)" != colnames(X), drop = FALSE]
    }

    # pass extracted data on to [cve.call()]
    dr <- cve.call(X, Y, method = method, max.dim = max.dim, ...)

    # overwrite `call` property from [cve.call()]
    dr$call <- match.call()
    return(dr)
}

#' @inherit cve title
#' @inherit cve description
#'
#' @param X Design predictor matrix.
#' @param Y \eqn{n}-dimensional vector of responses.
#' @param h bandwidth or function to estimate bandwidth, defaults to internaly
#'      estimated bandwidth.
#' @param nObs parameter for choosing bandwidth \code{h} using
#'   \code{\link{estimate.bandwidth}} (ignored if \code{h} is supplied).
#' @param method This character string specifies the method of fitting. The
#' options are
#' \itemize{
#'    \item \code{"mean"} method to estimate the mean subspace, see [1].
#'    \item \code{"central"} ensemble method to estimate the central subspace,
#'      see [2].
#'    \item \code{"weighted.mean"} variation of \code{"mean"} method with
#'      adaptive weighting of slices, see [1].
#'    \item \code{"weighted.central"} variation of \code{"central"} method with
#'      adaptive weighting of slices, see [2].
#' }
#' @param func_list a list of functions applied to \code{Y} used by ECVE
#'      (see [2]) for central subspace estimation. The default ensemble are
#'      indicator functions of the \eqn{[0, 10], (10, 20], ..., (90, 100]}
#'      percent response quantiles. (only relevant if \code{method} is
#'      \code{"central"} or \code{"weighted.central"}, ignored otherwise)
#' @param k Dimension of lower dimensional projection, if \code{k} is given
#'      only the specified dimension \code{B} matrix is estimated.
#' @param min.dim lower bounds for \code{k}, (ignored if \code{k} is supplied).
#' @param max.dim upper bounds for \code{k}, (ignored if \code{k} is supplied).
#' @param tau Initial step-size.
#' @param tol Tolerance for break condition.
#' @param max.iter maximum number of optimization steps.
#' @param attempts If \code{V.init} not supplied, the optimization is carried
#'      out \code{attempts} times with starting values drawn from the invariant
#'      measure on the Stiefel manifold (see \code{\link{rStiefel}}).
#' @param nr.proj The number of projection used for projective resampling for
#'      multivariate response \eqn{Y} (under active development, ignored for
#'      univariate response).
#' @param momentum number of \eqn{[0, 1)} giving the ration of momentum for
#'      eucledian gradient update with a momentum term. \code{momentum = 0}
#'      corresponds to normal gradient descend.
#' @param slack Positive scaling to allow small increases of the loss while
#'      optimizing, i.e. \code{slack = 0.1} allows the target function to
#'      increase up to \eqn{10 \%} in one optimization step.
#' @param gamma step-size reduction multiple. If gradient step with step size
#'      \code{tau} is not accepted \code{gamma * tau} is set to the next step
#'      size.
#' @param V.init Semi-orthogonal matrix of dimensions `(ncol(X), ncol(X) - k)
#'      used as starting value in the optimization. (If supplied,
#'      \code{attempts} is set to 0 and \code{k} to match dimension).
#' @param logger a logger function (only for advanced users, slows down the
#'      computation).
#'
#' @inherit cve return
#' @inherit cve references
#'
#' @examples
#' # create B for simulation (k = 1)
#' B <- rep(1, 5) / sqrt(5)
#'
#' set.seed(21)
#' # creat predictor data X ~ N(0, I_p)
#' X <- matrix(rnorm(500), 100, 5)
#' # simulate response variable
#' #     Y = f(B'X) + err
#' # with f(x1) = x1 and err ~ N(0, 0.25^2)
#' Y <- X %*% B + 0.25 * rnorm(100)
#'
#' # calculate cve with method 'simple' for k = 1
#' set.seed(21)
#' cve.obj.simple1 <- cve(Y ~ X, k = 1)
#'
#' # same as
#' set.seed(21)
#' cve.obj.simple2 <- cve.call(X, Y, k = 1)
#'
#' # extract estimated B's.
#' coef(cve.obj.simple1, k = 1)
#' coef(cve.obj.simple2, k = 1)
#' @export
cve.call <- function(X, Y,
    method = c("mean", "weighted.mean", "central", "weighted.central"),
    func_list = NULL, nObs = sqrt(nrow(X)), h = NULL,
    min.dim = 1L, max.dim = 10L, k = NULL,
    momentum = 0.0, tau = 1.0, tol = 1e-3,
    slack = 0.0, gamma = 0.5,
    V.init = NULL,
    max.iter = 50L, attempts = 10L, nr.proj = 1L,
    logger = NULL
) {
    # Determine method with partial matching (shortcuts: "Weight" -> "weighted")
    method <- match.arg(method)
    method_nr <- if(startsWith(method, "weighted")) 1L else 0L

    # Set default functions for ensamble methods (of indentity else)
    if (is.null(func_list)) {
        if (endsWith(method, "central")) {
            func_list <- list(
                function(Y) { q <- quantile(Y, 0.1); as.double(Y <= q) },
                function(Y) { q <- quantile(Y, c(0.1, 0.2)); as.double(q[1] < Y & Y <= q[2]) },
                function(Y) { q <- quantile(Y, c(0.2, 0.3)); as.double(q[1] < Y & Y <= q[2]) },
                function(Y) { q <- quantile(Y, c(0.3, 0.4)); as.double(q[1] < Y & Y <= q[2]) },
                function(Y) { q <- quantile(Y, c(0.4, 0.5)); as.double(q[1] < Y & Y <= q[2]) },
                function(Y) { q <- quantile(Y, c(0.5, 0.6)); as.double(q[1] < Y & Y <= q[2]) },
                function(Y) { q <- quantile(Y, c(0.6, 0.7)); as.double(q[1] < Y & Y <= q[2]) },
                function(Y) { q <- quantile(Y, c(0.7, 0.8)); as.double(q[1] < Y & Y <= q[2]) },
                function(Y) { q <- quantile(Y, c(0.8, 0.9)); as.double(q[1] < Y & Y <= q[2]) },
                function(Y) { q <- quantile(Y, 0.9); as.double(q < Y) }
            )
        } else {
            func_list <- list(function(Y) Y)
        }
    }

    # parameter checking
    if (!is.numeric(momentum) || length(momentum) > 1L) {
        stop("Momentum must be a number.")
    }
    if (!is.double(momentum)) {
        momentum <- as.double(momentum)
    }
    if (momentum < 0.0 || momentum >= 1.0) {
        stop("Momentum must be in [0, 1).")
    }

    if (!(is.matrix(X) && is.numeric(X))) {
        stop("Parameter 'X' should be a numeric matrix.")
    }
    if (!is.numeric(Y)) {
        stop("Parameter 'Y' must be numeric.")
    }
    if (!is.double(Y)) {
        storage.mode(Y) <- "double"
    }
    if (!is.matrix(Y)) {
        Y <- as.matrix(Y)
    }
    if (nrow(X) != nrow(Y)) {
        stop("Rows of 'X' and 'Y' elements are not compatible.")
    }
    if (ncol(X) < 2) {
        stop("'X' is one dimensional, no need for dimension reduction.")
    }

    if (!is.null(V.init)) {
        if (!is.matrix(V.init)) {
            stop("'V.init' must be a matrix.")
        }
        if (!all.equal(crossprod(V.init), diag(1, ncol(V.init)))) {
            stop("'V.init' must be semi-orthogonal.")
        }
        if (ncol(X) != nrow(V.init) || ncol(X) <= ncol(V.init)) {
            stop("Dimension missmatch of 'V.init' and 'X'")
        }
        min.dim <- max.dim <- ncol(X) - ncol(V.init)
        storage.mode(V.init) <- "double"
        attempts <- 0L
    } else if (missing(k) || is.null(k)) {
        min.dim <- as.integer(min.dim)
        max.dim <- as.integer(min(max.dim, ncol(X) - 1L))
    } else {
        min.dim <- max.dim <- as.integer(k)
    }
    if (min.dim > max.dim) {
        stop("'min.dim' bigger 'max.dim'.")
    }
    if (max.dim >= ncol(X)) {
        stop("'max.dim' (or 'k') must be smaller than 'ncol(X)'.")
    }

    if (missing(h) || is.null(h)) {
        estimate <- TRUE
    } else if (is.function(h)) {
        estimate <- TRUE
        estimate.bandwidth <- h
    } else if (is.numeric(h) && h > 0.0) {
        estimate <- FALSE
        h <- as.double(h)
    } else {
        stop("Bandwidth 'h' must be positive numeric.")
    }

    if (!is.numeric(tau) || length(tau) > 1L || tau <= 0.0) {
        stop("Initial step-width 'tau' must be positive number.")
    } else {
        tau <- as.double(tau)
    }
    if (!is.numeric(tol) || length(tol) > 1L || tol < 0.0) {
        stop("Break condition tolerance 'tol' must be not negative number.")
    } else {
        tol <- as.double(tol)
    }
    if (!is.numeric(slack) || length(slack) > 1L) {
        stop("Break condition slack 'slack' must be not negative number.")
    } else {
        slack <- as.double(slack)
    }
    if (!is.numeric(gamma) || length(gamma) > 1L || gamma <= 0.0 || gamma >= 1.0) {
        stop("Stepsize reduction 'gamma' must be between 0 and 1.")
    } else {
        gamma <- as.double(gamma)
    }

    if (!is.numeric(max.iter) || length(max.iter) > 1L) {
        stop("Parameter 'max.iter' must be positive integer.")
    } else if (!is.integer(max.iter)) {
        max.iter <- as.integer(max.iter)
    }
    if (max.iter < 1L) {
        stop("Parameter 'max.iter' must be at least 1L.")
    }

    if (!is.integer(nr.proj)) {
        nr.proj <- as.integer(nr.proj)
    }
    if (length(nr.proj) > 1 || nr.proj < 1) {
        stop("Parameter 'nr.proj' must be a single positive integer.")
    }

    if (is.null(V.init)) {
        if (!is.numeric(attempts) || length(attempts) > 1L) {
            stop("Parameter 'attempts' must be positive integer.")
        } else if (!is.integer(attempts)) {
            attempts <- as.integer(attempts)
        }
        if (attempts < 1L) {
            stop("Parameter 'attempts' must be at least 1L.")
        }
    }

    # Projective resampling of the multivariate `Y` as a `n x nr.proj` matrix.
    if (ncol(Y) > 1) {
        projections <- matrix(rnorm(ncol(Y) * nr.proj), nr.proj)
        projections <- projections / sqrt(rowSums(projections^2))
        Y <- Y %*% t(projections)
    }
    # Evaluate each function given (possible projected) `Y` and build a
    # `n x (|func_list| * nr.proj)` matrix.
    Fy <- vapply(func_list, do.call, Y, list(Y))
    dim(Fy) <- c(nrow(Fy), prod(dim(Fy)[-1]))

    # Convert numerical values to "double".
    storage.mode(X) <- storage.mode(Fy) <- "double"

    if (is.function(logger)) {
        loggerEnv <- environment(logger)
    } else {
        loggerEnv <- NULL
    }

    # Call specified method.
    method <- tolower(method)
    call <- match.call()
    dr <- list()
    dr$res <- list()
    for (k in min.dim:max.dim) {

        if (estimate) {
            h <- estimate.bandwidth(X, k, nObs)
        }

        dr.k <- .Call('cve_export', PACKAGE = 'CVarE',
                      X, Fy, k, h,
                      method_nr,
                      V.init,
                      momentum, tau, tol,
                      slack, gamma,
                      max.iter, attempts,
                      logger, loggerEnv)

        dr.k$B <- null(dr.k$V)
        dr.k$loss <- mean(dr.k$L)
        dr.k$h <- h
        dr.k$k <- k
        class(dr.k) <- "cve.k"
        dr$res[[as.character(k)]] <- dr.k
    }

    # augment result information
    dr$X <- X
    dr$Y <- Y
    dr$Fy <- Fy
    dr$method <- method
    dr$call <- call
    class(dr) <- "cve"
    return(dr)
}
