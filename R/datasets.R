#' Multivariate Normal Distribution.
#'
#' Random generation for the multivariate normal distribution.
#' \deqn{X \sim N_p(\mu, \Sigma)}{X ~ N_p(\mu, \Sigma)}
#'
#' @param n number of samples.
#' @param mu mean
#' @param sigma covariance matrix.
#'
#' @return a \eqn{n\times p}{n x p} matrix with samples in its rows.
#'
#' @examples
#' CVarE:::rmvnorm(20, sigma = matrix(c(2, 1, 1, 2), 2))
#' CVarE:::rmvnorm(20, mu = c(3, -1, 2))
#'
#' @keywords internal
rmvnorm <- function(n = 1, mu = rep(0, p), sigma = diag(p)) {
    if (!missing(sigma)) {
        p <- nrow(sigma)
    } else if (!missing(mu)) {
        mu <- matrix(mu, ncol = 1)
        p <- nrow(mu)
    } else {
        stop("At least one of 'mu' or 'sigma' must be supplied.")
    }

    return(rep(mu, each = n) + matrix(rnorm(n * p), n) %*% chol(sigma))
}

#' Multivariate t Distribution.
#'
#' Random generation from multivariate t distribution (student distribution).
#'
#' @param n number of samples.
#' @param mu mean
#' @param sigma a \eqn{k\times k}{k x k} positive definite matrix. If the degree
#' \eqn{\nu} if bigger than 2 the created covariance is
#' \deqn{var(x) = \Sigma\frac{\nu}{\nu - 2}}
#' for \eqn{\nu > 2}.
#' @param df degree of freedom \eqn{\nu}.
#'
#' @return a \eqn{n\times p}{n x p} matrix with samples in its rows.
#'
#' @examples
#' CVarE:::rmvt(20, c(0, 1), matrix(c(3, 1, 1, 2), 2), 3)
#' CVarE:::rmvt(20, sigma = matrix(c(2, 1, 1, 2), 2), df = 3)
#' CVarE:::rmvt(20, mu = c(3, -1, 2), df = 3)
#'
#' @keywords internal
rmvt <- function(n = 1, mu = rep(0, p), sigma = diag(p), df = Inf) {
    if (!missing(sigma)) {
        p <- nrow(sigma)
    } else if (!missing(mu)) {
        mu <- matrix(mu, ncol = 1)
        p <- nrow(mu)
    } else {
        stop("At least one of 'mu' or 'sigma' must be supplied.")
    }

    if (df == Inf) {
        Z <- 1
    } else {
        Z <- sqrt(df / rchisq(n, df))
    }

    return(rmvnorm(n, sigma = sigma) * Z + rep(mu, each = n))
}

#' Generalized Normal Distribution.
#'
#' Random generation for generalized Normal Distribution.
#'
#' @param n Number of generated samples.
#' @param mu mean.
#' @param alpha first shape parameter.
#' @param beta second shape parameter.
#'
#' @return numeric array of length \eqn{n}.
#'
#' @keywords internal
rgnorm <- function(n = 1, mu = 0, alpha = 1, beta = 1) {
    if (alpha <= 0 | beta <= 0) {
        stop("alpha and beta must be positive.")
    }
    lambda <- (1 / alpha)^beta
    scales <- qgamma(runif(n), shape = 1 / beta, scale = 1 / lambda)^(1 / beta)
    return(scales * ((-1)^rbinom(n, 1, 0.5)) + mu)
}

#' Laplace distribution
#'
#' Random generation for Laplace distribution.
#'
#' @param n Number of generated samples.
#' @param mu mean.
#' @param sd standard deviation.
#'
#' @return numeric array of length \eqn{n}.
#'
#' @keywords internal
rlaplace <- function(n = 1, mu = 0, sd = 1) {
    U <- runif(n, -0.5, 0.5)
    scale <- sd / sqrt(2)

    return(mu - scale * sign(U) * log(1 - 2 * abs(U)))
}

#' Generates test datasets.
#'
#' Provides sample datasets M1-M7 used in the paper Conditional variance
#' estimation for sufficient dimension reduction, Lukas Fertl, Efstathia Bura.
#' The general model is given by:
#' \deqn{Y = g(B'X) + \epsilon}
#'
#' @param name One of \code{"M1"}, \code{"M2"}, \code{"M3"}, \code{"M4",}
#' \code{"M5"}, \code{"M6"} or \code{"M7"}. Alternative just the dataset number
#' 1-7.
#' @param n number of samples.
#' @param p Dimension of random variable \eqn{X}.
#' @param sd standard diviation for error term \eqn{\epsilon}.
#' @param ... Additional parameters only for "M2" (namely \code{pmix} and
#' \code{lambda}), see: below.
#'
#' @return List with elements
#' \itemize{
#'      \item{X}{data, a \eqn{n\times p}{n x p} matrix.}
#'      \item{Y}{response.}
#'      \item{B}{the dim-reduction matrix}
#'      \item{name}{Name of the dataset (name parameter)}
#' }
#'
#' @section M1:
#' The predictors are distributed as
#' \eqn{X\sim N_p(0, \Sigma)}{X ~ N_p(0, \Sigma)} with
#' \eqn{\Sigma_{i, j} = 0.5^{|i - j|}}{\Sigma_ij = 0.5^|i - j|} for
#' \eqn{i, j = 1,..., p} for a subspace dimension of \eqn{k = 1} with a default
#' of \eqn{n = 100} data points. \eqn{p = 20},
#' \eqn{b_1 = (1,1,1,1,1,1,0,...,0)' / \sqrt{6}\in\mathcal{R}^p}{b_1 = (1,1,1,1,1,1,0,...,0)' / sqrt(6)}, and \eqn{Y} is
#' given as \deqn{Y = cos(b_1'X) + \epsilon} where \eqn{\epsilon} is
#' distributed as generalized normal distribution with location 0,
#' shape-parameter 0.5, and the scale-parameter is chosen such that
#' \eqn{Var(\epsilon) = 0.5}.
#' @section M2:
#' The predictors are distributed as \eqn{X \sim Z 1_p \lambda + N_p(0, I_p)}{X ~ Z 1_p \lambda + N_p(0, I_p)}. with
#' \eqn{Z \sim 2 Binom(p_{mix}) - 1\in\{-1, 1\}}{Z~2Binom(pmix)-1} where
#' \eqn{1_p} is the \eqn{p}-dimensional vector of one's, for a subspace
#' dimension of \eqn{k = 1} with a default of \eqn{n = 100} data points.
#' \eqn{p = 20}, \eqn{b_1 = (1,1,1,1,1,1,0,...,0)' / \sqrt{6}\in\mathcal{R}^p}{b_1 = (1,1,1,1,1,1,0,...,0)' / sqrt(6)},
#' and \eqn{Y} is \deqn{Y = cos(b_1'X) + 0.5\epsilon} where \eqn{\epsilon} is
#' standard normal.
#' Defaults for \code{pmix} is 0.3 and \code{lambda} defaults to 1.
#' @section M3:
#' The predictors are distributed as \eqn{X\sim N_p(0, I_p)}{X~N_p(0, I_p)}
#' for a subspace
#' dimension of \eqn{k = 1} with a default of \eqn{n = 100} data points.
#' \eqn{p = 20}, \eqn{b_1 = (1,1,1,1,1,1,0,...,0)' / \sqrt{6}\in\mathcal{R}^p}{b_1 = (1,1,1,1,1,1,0,...,0)' / sqrt(6)},
#' and \eqn{Y} is 
#' \deqn{Y = 2 log(|b_1'X| + 2) + 0.5\epsilon} where \eqn{\epsilon} is
#' standard normal.
#' @section M4:
#' The predictors are distributed as \eqn{X\sim N_p(0,\Sigma)}{X~N_p(0,\Sigma)}
#' with \eqn{\Sigma_{i, j} = 0.5^{|i - j|}}{\Sigma_ij = 0.5^|i - j|} for
#' \eqn{i, j = 1,..., p} for a subspace dimension of \eqn{k = 2} with a default
#' of \eqn{n = 100} data points. \eqn{p = 20},
#' \eqn{b_1 = (1,1,1,1,1,1,0,...,0)' / \sqrt{6}\in\mathcal{R}^p}{b_1 = (1,1,1,1,1,1,0,...,0)' / sqrt(6)},
#' \eqn{b_2 = (1,-1,1,-1,1,-1,0,...,0)' / \sqrt{6}\in\mathcal{R}^p}{b_2 = (1,-1,1,-1,1,-1,0,...,0)' / sqrt(6)}
#' and \eqn{Y} is given as \deqn{Y = \frac{b_1'X}{0.5 + (1.5 + b_2'X)^2} + 0.5\epsilon}{Y = (b_1'X) / (0.5 + (1.5 + b_2'X)^2) + 0.5\epsilon}
#' where \eqn{\epsilon} is standard normal.
#' @section M5:
#' The predictors are distributed as \eqn{X\sim U([0,1]^p)}{X~U([0, 1]^p)}
#' where \eqn{U([0, 1]^p)} is the uniform distribution with
#' independent components on the \eqn{p}-dimensional hypercube for a subspace
#' dimension of \eqn{k = 2} with a default of \eqn{n = 200} data points.
#' \eqn{p = 20},
#' \eqn{b_1 = (1,1,1,1,1,1,0,...,0)' / \sqrt{6}\in\mathcal{R}^p}{b_1 = (1,1,1,1,1,1,0,...,0)' / sqrt(6)},
#' \eqn{b_2 = (1,-1,1,-1,1,-1,0,...,0)' / \sqrt{6}\in\mathcal{R}^p}{b_2 = (1,-1,1,-1,1,-1,0,...,0)' / sqrt(6)}
#' and \eqn{Y} is given as \deqn{Y = cos(\pi b_1'X)(b_2'X + 1)^2 + 0.5\epsilon}
#' where \eqn{\epsilon} is standard normal.
#' @section M6:
#' The predictors are distributed as \eqn{X\sim N_p(0, I_p)}{X~N_p(0, I_p)}
#' for a subspace dimension of \eqn{k = 3} with a default of \eqn{n = 200} data
#' point. \eqn{p = 20, b_1 = e_1, b_2 = e_2}, and \eqn{b_3 = e_p}, where
#' \eqn{e_j} is the \eqn{j}-th unit vector in the \eqn{p}-dimensional space.
#' \eqn{Y} is given as \deqn{Y = (b_1'X)^2+(b_2'X)^2+(b_3'X)^2+0.5\epsilon}
#' where \eqn{\epsilon} is standard normal.
#' @section M7:
#' The predictors are distributed as \eqn{X\sim t_3(I_p)}{X~t_3(I_p)} where
#' \eqn{t_3(I_p)} is the standard multivariate t-distribution with 3 degrees of
#' freedom, for a subspace dimension of \eqn{k = 4} with a default of
#' \eqn{n = 200} data points.
#' \eqn{p = 20, b_1 = e_1, b_2 = e_2, b_3 = e_3}, and \eqn{b_4 = e_p}, where
#' \eqn{e_j} is the \eqn{j}-th unit vector in the \eqn{p}-dimensional space.
#' \eqn{Y} is given as \deqn{Y = (b_1'X)(b_2'X)^2+(b_3'X)(b_4'X)+0.5\epsilon}
#' where \eqn{\epsilon} is distributed as generalized normal distribution with
#' location 0, shape-parameter 1, and the scale-parameter is chosen such that
#' \eqn{Var(\epsilon) = 0.25}.
#'
#' @references
#'    Fertl, L. and Bura, E. (2021) "Conditional Variance
#'      Estimation for Sufficient Dimension Reduction"
#'      <arXiv:2102.08782>
#'
#' @import stats
#' @importFrom stats rnorm rbinom
#' @export
dataset <- function(name = "M1", n = NULL, p = 20, sd = 0.5, ...) {
    name <- toupper(name)
    if (nchar(name) == 1) { name <- paste0("M", name) }

    if (name == "M1") {
        if (missing(n)) { n <- 100 }
        # B ... `p x 1`
        B <- matrix(c(rep(1 / sqrt(6), 6), rep(0, p - 6)), ncol = 1)
        X <- rmvnorm(n, sigma = 0.5^abs(outer(1:p, 1:p, FUN = `-`)))
        beta <- 0.5
        Y <- cos(X %*% B) + rgnorm(n, 0,
            alpha = sqrt(sd^2 * gamma(1 / beta) / gamma(3 / beta)),
            beta = beta
        )
    } else if (name == "M2") {
        if (missing(n)) { n <- 100 }
        params <- list(...)
        pmix <- if (is.null(params$pmix)) { 0.3 } else { params$pmix }
        lambda <- if (is.null(params$lambda)) { 1 } else { params$lambda }
        # B ... `p x 1`
        B <- matrix(c(rep(1 / sqrt(6), 6), rep(0, p - 6)), ncol = 1)
        Z <- 2 * rbinom(n, 1, pmix) - 1
        X <- matrix(rep(lambda * Z, p) + rnorm(n * p), n)
        Y <- cos(X %*% B) + rnorm(n, 0, sd)
    } else if (name == "M3") {
        if (missing(n)) { n <- 100 }
        # B ... `p x 1`
        B <- matrix(c(rep(1 / sqrt(6), 6), rep(0, p - 6)), ncol = 1)
        X <- matrix(rnorm(n * p), n)
        Y <- 2 * log(2 + abs(X %*% B)) + rnorm(n, 0, sd)
    } else if (name == "M4") {
        if (missing(n)) { n <- 200 }
        # B ... `p x 2`
        B <- cbind(
            c(rep(1 / sqrt(6), 6), rep(0, p - 6)),
            c(rep(c(1, -1), 3) / sqrt(6), rep(0, p - 6))
        )
        X <- rmvnorm(n, sigma = 0.5^abs(outer(1:p, 1:p, FUN = `-`)))
        XB <- X %*% B
        Y <- (XB[, 1]) / (0.5 + (XB[, 2] + 1.5)^2) + rnorm(n, 0, sd)
    } else if (name == "M5") {
        if (missing(n)) { n <- 200 }
        # B ... `p x 2`
        B <- cbind(
            c(rep(1,        6), rep(0, p - 6)),
            c(rep(c(1, -1), 3), rep(0, p - 6))
        ) / sqrt(6)
        X <- matrix(runif(n * p), n)
        XB <- X %*% B
        Y <- cos(XB[, 1] * pi) * (XB[, 2] + 1)^2 + rnorm(n, 0, sd)
    } else if (name == "M6") {
        if (missing(n)) { n <- 200 }
        # B ... `p x 3`
        B <- diag(p)[, -(3:(p - 1))]
        X <- matrix(rnorm(n * p), n)
        Y <- rowSums((X %*% B)^2) + rnorm(n, 0, sd)
    } else if (name == "M7") {
        if (missing(n)) { n <- 400 }
        # B ... `p x 4`
        B <- diag(p)[, -(4:(p - 1))]
        # "R"andom "M"ulti"V"ariate "S"tudent
        X <- rmvt(n = n, sigma = diag(p), df = 3)
        XB <- X %*% B
        Y <- (XB[, 1]) * (XB[, 2])^2 + (XB[, 3]) * (XB[, 4])
        Y <- Y + rlaplace(n, 0, sd)
    } else {
        stop("Got unknown dataset name.")
    }

    return(list(X = X, Y = Y, B = B, name = name))
}
