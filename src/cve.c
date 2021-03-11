#include <R_ext/Utils.h> // for R_CheckUserInterrupt

#include "cve.h"

void cve(const mat *X, const mat *Fy, const double h,
         const unsigned int method,
         const double momentum,
         const double tau_init, const double tol_init,
         const double slack, const double gamma,
         const int maxIter, const int attempts,
         mat *V, mat *L,
         SEXP logger, SEXP loggerEnv) {

    int n = X->nrow, p = X->ncol, q = V->ncol;
    int attempt = 0, iter;
    double loss, loss_last, loss_best, err, tau;
    double tol = tol_init * sqrt((double)(2 * q));
    double agility = -2.0 * (1.0 - momentum) / (h * h);
    double sumK;
    double c = agility / (double)n;

    /* Create further intermediate or internal variables. */
    mat *lvecD_e  = (void*)0;
    mat *Fy_sq    = (void*)0;
    mat *XV       = (void*)0;
    mat *lvecD    = (void*)0;
    mat *D        = (void*)0;
    mat *lvecK    = (void*)0;
    mat *K        = (void*)0;
    mat *colSumsK = (void*)0;
    mat *rowSumsL = (void*)0;
    mat *W        = (void*)0;
    mat *y1       = (void*)0;
    mat *y2       = (void*)0;
    mat *S        = (void*)0;
    mat *tmp1     = (void*)0;
    mat *tmp2     = (void*)0;
    mat *G        = (void*)0;
    mat *A        = (void*)0;
    mat *V_tau    = (void*)0;
    mat *V_best   = (void*)0;
    mat *L_best   = (void*)0;

    /* Allocate appropriate amount of working memory. */
    int workLen = 2 * (p + 1) * p;
    if (workLen < n) {
        workLen = n;
    }
    double *workMem = (double*)R_alloc(workLen, sizeof(double));

    lvecD_e = rowDiffSquareSums(X, lvecD_e);
    Fy_sq = hadamard(1.0, Fy, Fy, 0.0, Fy_sq);

    do {
        /* (Re)set learning rate. */
        tau = tau_init;

        /* Check if start value for `V` was supplied. */
        if (attempts > 0) {
            /* Sample start value from Stiefel manifold. */
            V = rStiefel(p, q, V, workMem);
        }

        /* Embed X_i's in V space */
        XV = matrixprod(1.0, X, V, 0.0, XV);
        /* Compute embedded distances */
        lvecD = lincomb(1.0, lvecD_e, -1.0, rowDiffSquareSums(XV, lvecD));
        /* Apply kernel to distances. */
        lvecK = applyKernel(lvecD, h, gauss, lvecK);
        /* Transform lower vectors lvecD, lvecK into sym. matrices. */
        D = lvecToSym(lvecD, 0.0, D);
        K = lvecToSym(lvecK, 1.0, K);
        /* Compute column sums of kernel matrix K */
        colSumsK = colSums(K, colSumsK);
        /* Normalize K columns to obtain weight matrix W */
        W = colApply(K, '/', colSumsK, W);
        /* first and second order weighted responses */
        y1 = matrixprod(1.0, W, Fy,    0.0, y1);
        y2 = matrixprod(1.0, W, Fy_sq, 0.0, y2);
        /* Compute losses */
        L = hadamard(-1.0, y1, y1, 1.0, copy(y2, L));
        /* Compute initial loss */
        if (method == weighted) {
            colSumsK = elemApply(colSumsK, '-', 1.0, colSumsK);
            sumK = sum(colSumsK);
            if (L->ncol == 1) {
                loss_last = dot(L, '*', colSumsK) / sumK;
            } else {
                loss_last = dot(rowSums(L, rowSumsL), '*', colSumsK) / sumK;
            }
            c = agility / sumK;
            /* Calculate the scaling matrix S */
            S = laplace(adjacence(L, Fy, y1, D, K, gauss, S), workMem);
        } else { /* simple */
            loss_last = mean(L);
            /* Calculate the scaling matrix S */
            S = laplace(adjacence(L, Fy, y1, D, W, gauss, S), workMem);
        }
        /* Gradient */
        tmp1 = matrixprod(1.0, S, X, 0.0, tmp1);
        tmp2 = crossprod(1.0, X, tmp1, 0.0, tmp2);
        G = matrixprod(c, tmp2, V, 0.0, G);

        if (logger) {
            callLogger(logger, loggerEnv, attempt, /* iter <- 0L */ -1,
                       L, V, G, loss_last, /* err <- NA */ -1.0, tau);
        }

        /* Compute Skew-Symmetric matrix `A` used in Cayley transform.
         * `A <- tau * (G V^T - V G^T) + 0 * A`*/
        A = skew(tau, G, V, 0.0, A);

        for (iter = 0; iter < maxIter; ++iter) {
            /* Before next iteration, check if the User has requested an
             * interrupt (aka. ^C, or "Stop" button).
             * If interrupted the algorithm will be exited here and everything
             * will be discharted! */
            R_CheckUserInterrupt();

            /* Move `V` along the gradient direction. */
            V_tau = cayleyTransform(A, V, V_tau, workMem);

            /* Embed X_i's in V space */
            XV = matrixprod(1.0, X, V_tau, 0.0, XV);
            /* Compute embedded distances */
            lvecD = lincomb(1.0, lvecD_e, -1.0, rowDiffSquareSums(XV, lvecD));
            /* Apply kernel to distances. */
            lvecK = applyKernel(lvecD, h, gauss, lvecK);
            /* Transform lower vectors lvecD, lvecK into sym. matrices. */
            D = lvecToSym(lvecD, 0.0, D);
            K = lvecToSym(lvecK, 1.0, K);
            /* Compute column sums of kernel matrix K */
            colSumsK = colSums(K, colSumsK);
            /* Normalize K columns to obtain weight matrix W */
            W = colApply(K, '/', colSumsK, W);
            /* first and second order weighted responses */
            y1 = matrixprod(1.0, W, Fy,    0.0, y1);
            y2 = matrixprod(1.0, W, Fy_sq, 0.0, y2);
            /* Compute losses */
            L = hadamard(-1.0, y1, y1, 1.0, copy(y2, L));
            /* Compute loss */
            if (method == weighted) {
                colSumsK = elemApply(colSumsK, '-', 1.0, colSumsK);
                sumK = sum(colSumsK);
                if (L->ncol == 1) {
                    loss = dot(L, '*', colSumsK) / sumK;
                } else {
                    loss = dot(rowSums(L, rowSumsL), '*', colSumsK) / sumK;
                }
            } else { /* simple */
                loss = mean(L);
            }

            /* Check if step is appropriate, iff not reduce learning rate. */
            if ((loss - loss_last) > loss_last * slack) {
                tau *= gamma;
                iter -= 1;
                A = elemApply(A, '*', gamma, A); // scale A by gamma
                continue;
            } else {
                tau /= gamma;
            }

            /* Compute error, use workMem. */
            err = dist(V, V_tau);

            /* Shift next step to current step and store loss to last. */
            V = copy(V_tau, V);
            loss_last = loss;

            if (logger) {
                callLogger(logger, loggerEnv, attempt, iter,
                           L, V, G, loss, err, tau);
            }

            /* Check Break condition. */
            if (err < tol || iter + 1 >= maxIter) {
                break;
            }

            if (method == weighted) {
                /* Calculate the scaling matrix S */
                S = laplace(adjacence(L, Fy, y1, D, K, gauss, S), workMem);
                c = agility / sumK; // n removed previously
            } else { /* simple */
                /* Calculate the scaling matrix S */
                S = laplace(adjacence(L, Fy, y1, D, W, gauss, S), workMem);
            }

            /* Gradient */
            tmp1 = matrixprod(1.0, S, X, 0.0, tmp1);
            tmp2 = crossprod(1.0, X, tmp1, 0.0, tmp2);
            G = matrixprod(c, tmp2, V, momentum, G);

            /* Compute Skew-Symmetric matrix `A` used in Cayley transform.
             * `A <- tau * (G V^T - V G^T) + 0 * A`*/
            A = skew(tau, G, V, 0.0, A);
        }

        /* Check if current attempt improved previous ones */
        if (attempt == 0 || loss < loss_best) {
            loss_best = loss;
            V_best = copy(V, V_best);
            L_best = copy(L, L_best);
        }
    } while (++attempt < attempts);

    V = copy(V_best, V);
    L = copy(L_best, L);
}
