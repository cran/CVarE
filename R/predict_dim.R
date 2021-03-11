predict_dim_cv <- function(object) {
    # Get centered training data and dimensions
    X <- scale(object$X, center = TRUE, scale = FALSE)
    n <- nrow(object$X) # umber of training data samples
    Sigma <- (1 / n) * crossprod(X, X)
    eig <- eigen(Sigma)
    Sigma_root <- eig$vectors %*% tcrossprod(diag(sqrt(eig$values)), eig$vectors)
    X <- X %*% solve(Sigma_root)

    pred <- array(0, c(n, ncol(object$Fy), length(object$res)),
                  dimnames = list(NULL, NULL, names(object$res)))
    for (dr.k in object$res) {
        # get "name" of current dimension
        k <- as.character(dr.k$k)
        # Project dataset with current SDR basis
        X.proj <- X %*% dr.k$B

        for (i in 1:n) {
            model <- mda::mars(X.proj[-i, ], object$Fy[-i, ])
            pred[i, , k] <- predict(model, X.proj[i, , drop = FALSE])
        }
    }
    MSE <- apply((pred - as.numeric(object$Fy))^2, 3, mean)

    return(list(
        MSE = MSE,
        k = as.integer(names(which.min(MSE)))
    ))
}

predict_dim_elbow <- function(object) {
    if (ncol(object$Fy) > 1) # TODO: Implement or find better way
        stop("For multivariate or central models not supported yet.")

    # extract original data from object (cve result)
    X <- object$X
    Y <- object$Y
    # Get dimensions
    n <- nrow(X)
    p <- ncol(X)

    losses <- vector("double", length(object$res))
    names(losses) <- names(object$res)
    # Compute per sample losses with alternative bandwidth for each dimension.
    for (dr.k in object$res) {
        # extract dimension specific estimates and dimensions.
        k <- dr.k$k
        V <- dr.k$V
        # estimate bandwidth according alternative formula.
        h <- estimate.bandwidth(X, k, sqrt(n), version = 2L)
        # Projected `X`
        XQ <- X %*% (diag(1, p) - tcrossprod(V)) # X (I - V V')
        # Compute distances
        d2 <- tcrossprod(XQ) # XQ XQ'
        d1 <- matrix(diag(d2), n, n)
        D <- d1 - 2 * d2 + t(d1)
        # Apply kernel
        # Note: CVE uses for d = ||Q(X_i - X_j)|| the kernel exp(-d^4 / (2 h^2))
        K <- exp((-0.5 / h^2) * D^2)
        # sum columns
        colSumsK <- colSums(K)
        # compute weighted and square meighted reponses
        y1 <- (K %*% Y) / colSumsK
        y2 <- (K %*% Y^2) / colSumsK
        # total loss
        losses[[as.character(k)]] <- mean(y2 - y1^2)
    }

    return(list(
        losses = losses,
        k = as.integer(names(which.min(losses)))
    ))
}

predict_dim_wilcoxon <- function(object, p.value = 0.05) {
    if (ncol(object$Fy) > 1) # TODO: Implement or find better way
        stop("For multivariate or central models not supported yet.")

    # extract original data from object (cve result)
    X <- object$X
    Y <- object$Y
    # Get dimensions
    n <- nrow(X)
    p <- ncol(X)
    
    L <- matrix(NA, n, length(object$res))
    colnames(L) <- names(object$res)
    # Compute per sample losses with alternative bandwidth for each dimension.
    for (dr.k in object$res) {
        # extract dimension specific estimates and dimensions.
        k <- dr.k$k
        V <- dr.k$V
        # estimate bandwidth according alternative formula.
        h <- estimate.bandwidth(X, k, sqrt(n), version = 2L)
        # Projected `X`
        XQ <- X %*% (diag(1, p) - tcrossprod(V)) # X (I - V V')
        # Compute distances
        d2 <- tcrossprod(XQ) # XQ XQ'
        d1 <- matrix(diag(d2), n, n)
        D <- d1 - 2 * d2 + t(d1)
        # Apply kernel
        # Note: CVE uses for d = ||Q(X_i - X_j)|| the kernel exp(-d^4 / (2 h^2))
        K <- exp((-0.5 / h^2) * D^2)
        # sum columns
        colSumsK <- colSums(K)
        # compute weighted and square meighted reponses
        y1 <- (K %*% Y) / colSumsK
        y2 <- (K %*% Y^2) / colSumsK
        # element-wise L for dim. k
        L[, as.character(k)] <- y2 - y1^2
    }

    for (ind in seq_len(length(object$res) - 1L)) {
        p.test <- wilcox.test(L[, ind], L[, ind + 1L],
                              alternative = "less")$p.value
        if (p.test < p.value) {
            return(list(
                p.value = p.test,
                k = object$res[[ind]]$k
            ))
        }
    }

    return(list(
        p.value = NA,
        k = object$res[[length(object$res)]]$k
    ))
}

#' Estimate Dimension of the Sufficient Reduction.
#' 
#' This function estimates the dimension, i.e. the rank of \eqn{B}. The default
#' method \code{'CV'} performs leave-one-out (LOO) cross-validation using
#' \code{mars} as follows for \code{k = min.dim, ..., max.dim} a
#' cross-validation via \code{mars} is
#' performed on the dataset \eqn{(Y_i, B_k' X_i)_{i = 1, ..., n}} where
#' \eqn{B_k} is the \eqn{p \times k}{p x k} dimensional CVE estimate. The
#' estimated SDR dimension is the \eqn{k} where the
#' cross-validation mean squared error is minimal. The method \code{'elbow'}
#' estimates the dimension via \eqn{k = argmin_k L_n(V_{p - k})} where
#' \eqn{V_{p - k}} is the space that is orthogonal to the column space of the 
#' CVE estimate of \eqn{B_k}. Method \code{'wilcoxon'}  finds the minimum using
#' the Wilcoxon test.
#' 
#' @param object an object of class \code{"cve"}, usually, a result of a call to
#'          \code{\link{cve}} or \code{\link{cve.call}}.
#' @param method This parameter specifies which method is used in dimension
#'  estimation. It provides three options: \code{'CV'} (default),
#' \code{'elbow'} and \code{'wilcoxon'}.
#' @param ... ignored.
#'
#' @return A \code{list} with
#' \describe{
#'    \item{}{criterion for method and \code{k = min.dim, ..., max.dim}.}
#'    \item{k}{estimated dimension is the minimizer of the criterion.}
#' }
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
#' # Calculate cve for unknown k between min.dim and max.dim.
#' cve.obj.simple <- cve(y ~ x)
#'
#' predict_dim(cve.obj.simple)
#'
#' @export
predict_dim <- function(object, ..., method = "CV") {
    # Check if there are dimensions to select.
    if (length(object$res) == 1L) {
        return(list(
            message = "Only one dim. estimated.",
            k = as.integer(names(object$res))
        ))
    }

    # Determine method "fuzzy".
    methods <- c("cv", "elbow", "wilcoxon")
    names(methods) <- methods
    method <- methods[[tolower(method), exact = FALSE]]
    if (is.null(method)) {
        stop('Unable to determine method.')
    }

    if (method == "cv") {
        return(predict_dim_cv(object))
    } else if (method == "elbow") {
        return(predict_dim_elbow(object))
    } else if (method == "wilcoxon") {
        return(predict_dim_wilcoxon(object))
    } else {
        stop("Unable to determine method.")
    }
}
