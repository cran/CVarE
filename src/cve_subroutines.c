#include "cve.h"

/**
 * Applies the requested kernel element-wise.
 *
 * @param A matrix to apply the kernel to.
 * @param h bandwidth parameter.
 * @param kernel the kernel to be used.
 * @param B (in/out) matrix `A` with element-wise applied kernel.
 *
 * @returns ether the passed `B`, or a new created matrix if `B` was NULL.
 */
mat* applyKernel(const mat* A, const double h, kernel kernel, mat* B) {
    int i, nn = A->nrow * A->ncol;
    int nn4 = 4 * (nn / 4);
    double scale;
    const double * restrict a;
    double * restrict b;

    if (!B) {
        B = matrix(A->nrow, A->ncol);
    } else if (nn != B->nrow * B->ncol) {
        // TODO: error handling!
    }

    a = A->elem;
    b = B->elem;
    switch (kernel) {
        case gauss:
            scale = -0.5 / (h * h);
            for (i = 0; i < nn4; i += 4) {
                b[i]     = exp(scale * a[i]     * a[i]);
                b[i + 1] = exp(scale * a[i + 1] * a[i + 1]);
                b[i + 2] = exp(scale * a[i + 2] * a[i + 2]);
                b[i + 3] = exp(scale * a[i + 3] * a[i + 3]);
            }
            for (; i < nn; ++i) {
                b[i] = exp(scale * a[i] * a[i]);
            }
            break;
        default:
            // TODO: error handling!
            break;
    }

    return B;
}

/**
 * Scaling matrix `S` defined as
 *      s_{i j} = (L_i - (Y_j - y1_i)^2) * d_{i j} * w_{i j}
 *
 * Mapping of vectors `L`, `y1` and `Y` combinations.
 *
 *                                       Y[j]
 *                                  ------ j ----->
 *                  s[0]    s[n]  . . .  s[jn]  . . .  s[(j-1)n]
 *                  s[1]   s[n+1] . . . s[jn+1] . . . s[(j-1)n+1]
 *               |    .       .   .        .    .          .
 *               |    .       .     .      .      .        .
 *                    .       .       .    .        .      .
 *  L[i], y1[i]  i  s[i]   s[n+i] . . . s[jn+i] . . . s[(j-1)n+i]
 *                    .       .   .        .    .          .
 *               |    .       .     .      .      .        .
 *               v    .       .       .    .        .      .
 *                 s[n-1] s[2n-1] . . .  s[n-1] . . .   s[nn-1]
 *
 * @param L per sample loss vector of (length `n`).
 * @param Y responses (length `n`).
 * @param y1 weighted responses (length `n`).
 * @param D distance matrix (dim. `n x n`).
 * @param W weight matrix (dim. `n x n`).
 * @param kernel the kernel to be used.
 *
 * @returns passed matrix `S` if not NULL, then a new `n x n` matrix is created.
 *
 * @example Basically equivalent to the following R function.
 * r_LS <- function(L, Y, y1, D, W) {
 *     # get dimension
 *     n <- length(L)
 *     # Indices
 *     i <- rep(1:n, n)
 *     j <- rep(1:n, each = n)
 *     # Compute S
 *     matrix(L[i] - (Y[j] - y1[i])^2, n) * W * D
 * }
 * 
 * @details mapping for indices used in blocked implementation.
 *
 * n ..... Dimensions of `S`, a `n x n` matrix.
 * B ..... BLOCK_SIZE
 * rB .... block reminder, aka rB = B % n.
 *
 *                                    Y[j]
 *                            0  -----  j  ----->  n
 *                          0 +--------------------+
 *                            | '                  |
 *                          ' | k      B x n       |
 *                          ' | v                  |
 *                          ' +--------------------+
 *            L[i], y1[i]   i | '                  |
 *                          ' | k      B x n       |
 *                          ' | v                  |
 *                          v +--------------------+
 *                            | k     rB x n       |
 *                          n +--------------------+
 */
// TODO: fix: cache misses in Y?!
mat* adjacence(const mat *mat_L, const mat *mat_Fy, const mat *mat_y1,
               const mat *mat_D, const mat *mat_W, kernel kernel,
               mat *mat_S) {
    int i, j, k, l, n = mat_L->nrow, m = mat_L->ncol;
    int block_size, block_batch_size;
    int max_size = 64 < n ? 64 : n; // Block Size set to 64

    double Y_j, t0, t1, t2, t3; // internal temp. values.
    double *Y, *L, *y1;
    double *D, *W, *S;

    if (!mat_S) {
        mat_S = zero(n, n);
    } else {
        memset(mat_S->elem, 0, n * n * sizeof(double));
    }

    for (l = 0; l < m; ++l) {
        Y  = mat_Fy->elem + l * n;
        L  = mat_L->elem + l * n;
        y1 = mat_y1->elem + l * n;
        for (i = 0; i < n; i += block_size) {
            /* Set blocks (left upper corner) */
            S = mat_S->elem + i;
            D = mat_D->elem + i;
            W = mat_W->elem + i;
            /* determine block size */
            block_size = n - i;
            if (block_size > max_size) {
                block_size = max_size;
            }
            block_batch_size = 4 * (block_size / 4);
            /* for each column in the block */
            for (j = 0; j < n; ++j, S += n, D += n, W += n) {
                Y_j = Y[j];
                /* iterate over block rows */
                for (k = 0; k < block_batch_size; k += 4) {
                    t0 = Y_j - y1[k];
                    t1 = Y_j - y1[k + 1];
                    t2 = Y_j - y1[k + 2];
                    t3 = Y_j - y1[k + 3];
                    S[k]     += (L[k]     - (t0 * t0)) * D[k]     * W[k];
                    S[k + 1] += (L[k + 1] - (t1 * t1)) * D[k + 1] * W[k + 1];
                    S[k + 2] += (L[k + 2] - (t2 * t2)) * D[k + 2] * W[k + 2];
                    S[k + 3] += (L[k + 3] - (t3 * t3)) * D[k + 3] * W[k + 3];
                }
                for (; k < block_size; ++k) {
                    t0 = Y_j - y1[k];
                    S[k] += (L[k] - (t0 * t0)) * D[k] * W[k];
                }
            }
            L  += block_size;
            y1 += block_size;
        }
    }

    if (m > 1) {
        S = mat_S->elem;
        for (i = 0; i < n * n; ++i) {
            S[i] /= m;
        }
    }

    return mat_S;
}
