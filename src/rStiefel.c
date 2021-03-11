#include "cve.h"

/**
 * Draws a sample from invariant measure on the Stiefel manifold \eqn{S(p, q)}.
 *
 * @param p row dimension
 * @param q column dimension
 * @param V (in/out) matrix of dimensions `p x q` or NULL.
 * @param workMem work space array of length greater-equal than `2pq + q`.
 *
 * @return Passed matrix `V` or new created if `V` is NULL.
 *
 * @example Performs the same operation as the following `R` code:
 *      V <- qr.Q(qr(matrix(rnorm(p * q, 0, 1), p, q)))
 * 
 * @details ATTENTION: The length of workMem must be at least `2pq + q`.
 */
mat* rStiefel(const int p, const int q, mat *V, double *workMem) {
    int i, j, info, workLen = 2 * p * q + q;
    int pq = p * q;
    double *v;

    if (!V) {
        V = matrix(p, q);
    } else if (V->nrow != p || V->ncol != q) {
        // TODO: error handling!
    }
    v = V->elem;

    GetRNGstate();
    for (i = 0; i < pq; ++i) {
        workMem[i] = norm_rand();
    }
    PutRNGstate();

    double *tau = workMem + pq;
    workLen -= pq + q;

    F77_CALL(dgeqrf)(&p, &q, workMem, &p, tau,
                     workMem + pq + q, &workLen, &info);

    for (j = 0; j < q; ++j) {
        for (i = 0; i < p; ++i) {
            v[p * j + i] = i == j ? 1.0 : 0.0;
        }
    }

    F77_NAME(dormqr)("L", "N", &p, &q, &q, workMem, &p, tau, V->elem, &p,
                     workMem + pq + q, &workLen, &info);

    return V;
}
