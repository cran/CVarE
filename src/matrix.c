#include "cve.h"

/**
 * Creates a matrix.
 *
 * @param nrow number of rows.
 * @param ncol number of columns.
 * 
 * @returns an uninitialized `nrow x ncol` matrix. 
 *
 * @details matrix elements are NOT initialized, its content is "random".
 */
mat* matrix(const int nrow, const int ncol) {
    mat* M = (mat*)R_alloc(1, sizeof(mat));
    M->nrow = nrow;
    M->ncol = ncol;
    M->origin = (void*)0;
    M->elem = (double*)R_alloc(nrow * ncol, sizeof(double));
    return M;
}

/**
 * Creates a zero matrix.
 *
 * @param nrow number of rows.
 * @param ncol number of columns.
 * 
 * @returns a `nrow x ncol` matrix with all elements equal 0.
 */
mat* zero(const int nrow, const int ncol) {
    mat* M = matrix(nrow, ncol);
    /* Set all emements zero */
    memset(M->elem, 0, nrow * ncol * sizeof(double));
    return M;
}

/**
 * Copies elements form `src` to `dest` matrix.
 *
 * @param dest (in/out) destination matrix.
 * @param src source matrix.
 *
 * @returns passed dest matrix.
 */
mat* copy(mat *src, mat *dest) {
    if (!dest) {
        dest = matrix(src->nrow, src->ncol);
    } else if (src->nrow != dest->nrow || src->ncol != dest->ncol) {
        // TODO: error handling!
    }
    memcpy(dest->elem, src->elem, dest->nrow * dest->ncol * sizeof(double));

    return dest;
}

/**
 * Sum of all elements.
 *
 * @param A matrix.
 *
 * @returns the sum of elements of `A`.
 */
double sum(const mat* A) {
    int i, nn = A->nrow * A->ncol;
    int nn4 = 4 * (nn / 4);
    double *a = A->elem;
    double s = 0.0;

    for (i = 0; i < nn4; i += 4) {
        s += a[i] + a[i + 1] + a[i + 2] + a[i + 3];
    }
    for (; i < nn; ++i) {
        s += a[i];
    }

    return s;
}

/**
 * Mean of all elements.
 * 
 * @param A matrix.
 * 
 * @returns the mean of elements of `A`.
 */
double mean(const mat* A) {
    int i, nn = A->nrow * A->ncol;
    int nn4 = 4 * (nn / 4);
    double *a = A->elem;
    double s = 0.0;

    for (i = 0; i < nn4; i += 4) {
        s += a[i] + a[i + 1] + a[i + 2] + a[i + 3];
    }
    for (; i < nn; ++i) {
        s += a[i];
    }

    return s / (double)nn;
}

#define CVE_DOT_ALG(op)                   \
    for (i = 0; i < nn4; i += 4) {        \
        s += ((a[i])     op (b[i]))       \
           + ((a[i + 1]) op (b[i + 1]))   \
           + ((a[i + 2]) op (b[i + 2]))   \
           + ((a[i + 3]) op (b[i + 3]));  \
    }                                     \
    for (; i < nn; ++i) {                 \
        s += (a[i]) op (b[i]);            \
    }

/**
 * Element-wise sum of `A` and `op(B)`.
 * 
 *      sum(A + B), for op = '+',
 *      sum(A * B), for op = '*',
 *      sum(A - B), for op = '-',
 *      sum(A / B), for op = '/'.
 * 
 * @param A matrix.
 * @param op one of '+', '-', '*' or '/'.
 * @param B matrix of the same size as `A`.
 * 
 * @returns sum of element-wise products.
 */
double dot(const mat *A, const char op, const mat *B) {
    int i, nn = A->nrow * A->ncol;
    int nn4 = 4 * (nn / 4);
    double *a = A->elem;
    double *b = B->elem;
    double s = 0.0;

    if (B->nrow * B->ncol != nn) {
        // TODO: error handling!
    }

    switch (op) {
        case '+':
            CVE_DOT_ALG(+)
            break;
        case '-':
            CVE_DOT_ALG(-)
            break;
        case '*':
            CVE_DOT_ALG(*)
            break;
        case '/':
            CVE_DOT_ALG(/)
            break;
        default:
            // TODO: error handling!
            break;
    }

    return s;
}

/**
 * Computes the sub-space distances according
 * 
 *      dist <- norm(A %*% t(A) - B %*% t(B), 'F')
 * 
 * @param A semi-orthogonal matrix of dim. `n x m` such that `n < m`.
 * @param B semi-orthogonal matrix of same dimensions as `A`.
 * 
 * @details The semi-orthogonality will NOT be validated!
 * 
 * @returns dist as describet above.
 */
// TODO: optimize!
double dist(const mat *A, const mat *B) {
    int i, j, k, n = A->nrow, m = A->ncol;
    double tmp1, tmp2, sum;
    double *a = A->elem, *b = B->elem;

    if (A->nrow != B->nrow || A->ncol != B->ncol) {
        // TODO: error handling!
    }

    sum = 0.0;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            tmp1 = tmp2 = 0.0;
            for (k = 0; k < m; ++k) {
                tmp1 += a[k * n + i] * a[k * n + j];
                tmp2 += b[k * n + i] * b[k * n + j];
            }
            sum += (tmp1 - tmp2) * (tmp1 - tmp2);
        }
    }

    return sqrt(sum);
}


/**
 * Sum of squared elements.
 *
 * @param A matrix.
 *
 * @returns the sum of squared elements of `A`.
 */
double squareSum(const mat* A) {
    int i, nn = A->nrow * A->ncol;
    int nn4 = 4 * (nn / 4);
    double *a = A->elem;
    double s = 0.0;

    for (i = 0; i < nn4; i += 4) {
        s += (a[i]     * a[i])
           + (a[i + 1] * a[i + 1])
           + (a[i + 2] * a[i + 2])
           + (a[i + 3] * a[i + 3]);
    }
    for (; i < nn; ++i) {
        s += a[i] * a[i];
    }

    return s;
}


/**
 * Row sums of a matrix `A`.
 *
 * @param A matrix of size `n x m`.
 * @param sums vector of length `n`.
 *
 * @returns sums if not NULL, else a new creates vector of length `n`.
 */
mat* rowSums(const mat *A, mat *sums) {
    int i, j, block_size, block_size_i,
        nrow = A->nrow, ncol = A->ncol;
    const double *a;
    const double *A_block = A->elem;
    const double *A_end = A->elem + nrow * ncol;
    double *sum;

    if (!sums) {
        sums = zero(nrow, 1);
    }
    sum = sums->elem;

    if (nrow > CVE_MEM_CHUNK_SIZE) {
        block_size = CVE_MEM_CHUNK_SIZE;
    } else {
        block_size = nrow;
    }

    // Iterate `(block_size_i, ncol)` submatrix blocks.
    for (i = 0; i < nrow; i += block_size_i) {
        // Reset `A` to new block beginning.
        a = A_block;
        // Take block size of eveything left and reduce to max size.
        block_size_i = nrow - i;
        if (block_size_i > block_size) {
            block_size_i = block_size;
        }
        // Compute first blocks column,
        for (j = 0; j < block_size_i; ++j) {
            sum[j] = a[j];
        }
        // and sum the following columns to the first one.
        for (a += nrow; a < A_end; a += nrow) {
            for (j = 0; j < block_size_i; ++j) {
                sum[j] += a[j];
            }
        }
        // Step one block forth.
        A_block += block_size_i;
        sum += block_size_i;
    }

    return sums;
}

/**
 * Column sums of a matrix `A`.
 *
 * @param A matrix of size `n x m`.
 * @param sums (out) vector of length `m`.
 *
 * @returns sums if not NULL, else a new creates vector of length `m`.
 * 
 * @details As a vector this the result `sums` must have matrix dim. `m x 1`.
 * Meaning that the result is a row vector.
 */
mat* colSums(const mat *A, mat *sums) {
    int i, j, nrow = A->nrow, ncol = A->ncol;
    int nrowb = 4 * (nrow / 4); // 4 * floor(nrow / 4)
    double *a = A->elem;
    double sum;

    if (!sums) {
        sums = matrix(ncol, 1);
    }

    for (j = 0; j < ncol; ++j) {
        sum = 0.0;
        for (i = 0; i < nrowb; i += 4) {
            sum += a[i]
                 + a[i + 1]
                 + a[i + 2]
                 + a[i + 3];
        }
        for (; i < nrow; ++i) {
            sum += a[i];
        }
        sums->elem[j] = sum;
        a += nrow;
    }

    return sums;
}

/**
 * Sum of squared row differences.
 *
 *      d_k = \sum_{l} (X_{i l} - X_{j l})^2
 *
 * For a input matrix `X` of dimension `n x p` the result is a vector of length
 * `n * (n - 1) / 2` (aka number of all pairs of `n` elements). The relation
 * between `i, j` and `k` is:
 * 
 *      k = j n + i - j (j + 1) / 2, for 0 <= j < i < n
 *
 * @param X input matrix of dim. `n x p`.
 * @param lvecD (out) output vector of length `n * (n - 1) / 2`.
 *
 * @returns lvecD if not NULL, otherwise a new created vector.
 *
 * @details this computes the packed vectorized lower triangular part of the
 * distance matrix of `X`, aka lvec(D).
 */
mat* rowDiffSquareSums(const mat *X, mat *lvecD) {
    int i, j, k, l, nrow = X->nrow, ncol = X->ncol;
    const double *x;
    double *d, diff;

    if (!lvecD) {
        lvecD = zero(nrow * (nrow - 1) / 2, 1);
    } else if (nrow * (nrow - 1) == 2 * lvecD->nrow && lvecD->ncol == 1) {
        memset(lvecD->elem, 0, (nrow * (nrow - 1) / 2) * sizeof(double));
    } else {
        // TODO: error handling!
    }

    d = lvecD->elem;
    for (l = 0; l < ncol; ++l) {
        x = X->elem + l * nrow;
        for (j = k = 0; j < nrow; ++j) {
            for (i = j + 1; i < nrow; ++i, ++k) {
                diff = x[i] - x[j];
                d[k] += diff * diff;
            }
        }
    }

    return lvecD;
}


#define COLAPPLY_ALG(op)                                                      \
    for (i = 0; i < nrow; i += block_size) {                                  \
        /* Set a, c to block's left upper corner and b to block height */     \
        a = A->elem + i;                                                      \
        b = B->elem + i;                                                      \
        c = C->elem + i;                                                      \
        /* Get block size as everyting left, then restrict size. */           \
        block_size = nrow - i;                                                \
        if (block_size > max_size) {                                          \
            block_size = max_size;                                            \
        }                                                                     \
        /* remove reminders by 4 from block_size */                           \
        block_batch_size = 4 * (block_size / 4);                              \
        /* Stay on i-index "height" and move "right" in j-direction */        \
        for (j = 0; j < ncol; ++j, a += nrow, c += nrow) {                    \
            /*  Iterate over the block column */                              \
            for (k = 0; k < block_batch_size; k += 4) {                       \
                c[k]     = a[k]     op b[k];                                  \
                c[k + 1] = a[k + 1] op b[k + 1];                              \
                c[k + 2] = a[k + 2] op b[k + 2];                              \
                c[k + 3] = a[k + 3] op b[k + 3];                              \
            }                                                                 \
            for (; k < block_size; ++k) {                                     \
                c[k] = a[k] op b[k];                                          \
            }                                                                 \
        }                                                                     \
    }

/* C[, j] = A[, j] op v for each j = 1 to ncol with op as one of +, -, *, / */

/**
 * Elment-wise operation for each row to vector.
 *      C <- A op B
 * using recycling (column-major order) with `op` ether '+', '-', '*' or '/'.
 *
 * @param A matrix of size `n x m`.
 * @param op type of operation, ether '+', '-', '*' or '/'.
 * @param B vector or length `m` to be sweeped over rows.
 * @param C (in/out) result after row sweep.
 *
 * @returns Passed (overwritten) `C` or new created if `C` is NULL.
 */
mat* colApply(mat *A, const char op, mat *B, mat *C) {
    int i, j, k, nrow = A->nrow, ncol = A->ncol;
    int max_size = CVE_MEM_CHUNK_SMALL < nrow ? CVE_MEM_CHUNK_SMALL : nrow;
    int block_size, block_batch_size;
    const double * restrict a, * restrict b;
    double * restrict c;

    /* Check if B length is the same as count A's rows. */
    if (nrow != B->nrow) {
        // TODO: error handling!
    }
    /* Check for existance or shape of C */
    if (!C) {
        C = matrix(nrow, ncol);
    } else if (C->nrow != nrow || C->ncol != ncol) {
        // TODO: error handling!
    }

    switch (op) {
        case '+':
            COLAPPLY_ALG(+)
            break;
        case '-':
            COLAPPLY_ALG(-)
            break;
        case '*':
            COLAPPLY_ALG(*)
            break;
        case '/':
            // TODO: check division by zero
            COLAPPLY_ALG(/)
            break;
        default:
            // TODO: error handling:
            break;
    }

    return C;
}
#undef COLAPPLY_ALG

/**
 * Element-wise applies `op` with scalar.
 * 
 *      B <- A + scalar,    for op = '+',
 *      B <- A - scalar,    for op = '-',
 *      B <- A * scalar,    for op = '*',
 *      B <- A / scalar,    for op = '/'.
 * 
 * @param A matrix.
 * @param op one of '+', '-', '*' or '/'.
 * @param scalar scalar as right hand side of operation.
 * @param B matrix of the same size as `A`. May even `A`.
 */
mat* elemApply(mat *A, const char op, const double scalar, mat *B) {
    int i, n = A->nrow * A->ncol;
    double *a, *b;

    if (!B) {
        B = matrix(A->nrow, A->ncol);
    }

    if (n != B->nrow * B->ncol || (op == '/' && scalar == 0.0)) {
        // TODO: error handling!
    }

    a = A->elem;
    b = B->elem;
    switch (op) {
        case '+':
            for (i = 0; i < n; ++i) {
                b[i] = a[i] + scalar;
            }
            break;
        case '-':
            for (i = 0; i < n; ++i) {
                b[i] = a[i] - scalar;
            }
            break;
        case '*':
            for (i = 0; i < n; ++i) {
                b[i] = a[i] * scalar;
            }
            break;
        case '/':
            for (i = 0; i < n; ++i) {
                b[i] = a[i] / scalar;
            }
            break;
        default:
            break;
    }

    return B;
}

/**
 * Linear combination of A with B,
 *      B <- alpha A + beta B
 *
 * @param alpha scaling factor for `A`.
 * @param A (in) first matrix.
 * @param beta scaling factor for `B`.
 * @param B (in/out) second matrix of the same length as `A`.
 *
 * @returns passed vector `B`.
 *
 * @details It is NOT allowed to pass NULL.
 * Also the number of elements must match, but it is legal to pass matrices of
 * different shape, the result may be unwanted but leagel.
 */
mat* lincomb(const double alpha, const mat *A, const double beta, mat *B) {
    int i, nn = A->nrow * A->ncol;
    int nnb = 4 * (nn / 4);
    double *restrict a = A->elem, *restrict b = B->elem;

    if (nn != B->nrow * B->ncol) {
        // TODO: error handling!
    }

    for (i = 0; i < nnb; i += 4) {
        b[i]     = alpha * a[i]     + beta * b[i];
        b[i + 1] = alpha * a[i + 1] + beta * b[i + 1];
        b[i + 2] = alpha * a[i + 2] + beta * b[i + 2];
        b[i + 3] = alpha * a[i + 3] + beta * b[i + 3];
    }
    for (; i < nn; ++i) {
        b[i] = alpha * a[i] + beta * b[i];
    }

    return B;
}

/**
 * Matrix matrix product.
 *      C <- alpha * A %*% B + beta * C
 *
 * @param alpha scaling factor for product of `A` with `B`.
 * @param A (in) matrix of dimension `n x k`.
 * @param B (in) matrix of dimension `k x m`.
 * @param beta scaling factor for original values in target matrix.
 * @param C (in/out) matrix of dimension `n x m` or NULL.
 *
 * @details if `C` is NULL, a new matrix of dimension `n x m` is created `beta`
 *      is set to 0.
 * 
 * @returns given `C` or a new created matrix if `C` was NULL or NULL on error.
 */
mat* matrixprod(const double alpha, const mat *A, const mat *B,
                double beta, mat* C) {
    /* Check dimensions. */
    if (A->ncol != B->nrow) {
        // TODO: error handling!
    }
    /* Check if `C` is given. */
    if (!C) {
        /* with `beta` zero, DGEMM ensures to overwrite `C`. */
        beta = 0.0;
        /* Create uninitialized matrix `C`. */
        C = matrix(A->nrow, B->ncol);
    } else {
        /* Check target dimensions. */
        if (C->nrow != A->nrow || C->ncol != B->ncol) {
            // TODO: error handling!
        }
    }

    /* DGEMM ... Dence GEneralized Matrix Matrix level 3 BLAS routine */
    F77_NAME(dgemm)("N", "N", &(A->nrow), &(B->ncol), &(A->ncol),
                    &alpha, A->elem, &(A->nrow), B->elem, &(B->nrow),
                    &beta, C->elem, &(C->nrow));

    return C;
}



/**
 * Cross-product.
 *      C <- alpha * t(A) %*% B + beta * C
 *      C <- alpha * crossprod(A, B) + beta * C
 *
 * @param alpha scaling factor for product of `A` with `B`.
 * @param A (in) matrix of dimension `k x n`.
 * @param B (in) matrix of dimension `k x m`.
 * @param beta scaling factor for original values in target matrix.
 * @param C (in/out) matrix of dimension `n x m` or NULL.
 *
 * @details if `C` is NULL, a new matrix of dimension `n x m` is created `beta`
 *      is set to 0.
 * 
 * @returns given `C` or a new created matrix if `C` was NULL or NULL on error.
 */
mat* crossprod(const double alpha, const mat *A, const mat *B,
                double beta, mat* C) {
    /* Check dimensions. */
    if (A->nrow != B->nrow) {
        // TODO: error handling!
    }
    /* Check if `C` is given. */
    if (!C) {
        /* with `beta` zero, DGEMM ensures to overwrite `C`. */
        beta = 0.0;
        /* Create uninitialized matrix `C`. */
        C = matrix(A->ncol, B->ncol);
    } else {
        /* Check target dimensions. */
        if (C->nrow != A->ncol || C->ncol != B->ncol) {
            // TODO: error handling!
        }
    }

    /* DGEMM ... Dence GEneralized Matrix Matrix level 3 BLAS routine */
    F77_NAME(dgemm)("T", "N", &(A->ncol), &(B->ncol), &(A->nrow),
                    &alpha, A->elem, &(A->nrow), B->elem, &(B->nrow),
                    &beta, C->elem, &(C->nrow));

    return C;
}

/**
 * Hadamard product (Schur product or element-wise product).
 *
 *      C <- alpha * A * B + beta * C
 *
 * @param alpha scaling factor for `A o B`.
 * @param A matrix.
 * @param B matrix of same dimensions as `A`.
 * @param beta scaling factor for `C`.
 * @param C (in/out) matrix of same dimensions as `A` or NULL.
 *
 * @returns Passed C, or if C is NULL a new created matrix.
 */
mat* hadamard(const double alpha, const mat* A, const mat* B,
              double beta, mat* C) {
    int i, nn = A->nrow * A->ncol;
    double *a, *b, *c;

    if (A->nrow != B->nrow || A->ncol != B->ncol) {
        // TODO: error handling!
    }

    if (!C) {
        beta = 0.0;
        C = zero(A->nrow, A->ncol);
    } else if (A->nrow != C->nrow || A->ncol != C->ncol) {
        // TOD0: error handling!
    }

    a = A->elem;
    b = B->elem;
    c = C->elem;
    if (alpha == 1.0 && beta == 0.0) {
        for (i = 0; i < nn; ++i) {
            c[i] = a[i] * b[i];
        }
    } else {
        for (i = 0; i < nn; ++i) {
            c[i] = alpha * a[i] * b[i] + beta * c[i];
        }
    }

    return C;
}

/**
 * Skew-Symmetric rank-2 matrixprod:
 *
 *      C <- alpha * (A %*% t(B) - B %*% t(A)) + beta * C
 *
 * @param alpha first scaling factor.
 * @param A (in) matrix of dimension `n x k`.
 * @param B (in) matrix of dimension `n x k`.
 * @param beta scaling factor for original values in target matrix.
 * @param C (in/out) matrix of dimension `n x n` or NULL.
 *
 * @details if `C` is NULL, a new matrix of dimension `n x n` is created, `beta`
 *      is set to 0.
 * 
 * @returns given `C` or a new created matrix if `C` was NULL or NULL on error.
 */
mat* skew(double alpha, mat *A, mat *B, double beta, mat *C) {

    if (A->nrow != B->nrow || A->ncol != B->ncol) {
        // TODO: error handling!
    }

    if (!C) {
        beta = 0.0;
        C = matrix(A->nrow, A->nrow);
    } else if (C->nrow != A->nrow || C->ncol != A->nrow) {
        // TODO: error handling!
    }

    F77_NAME(dgemm)("N", "T",
                    &(C->nrow), &(C->ncol), &(A->ncol),
                    &alpha, A->elem, &(A->nrow), B->elem, &(B->nrow),
                    &beta, C->elem, &(C->nrow));

    alpha *= -1.0;
    beta = 1.0;
    F77_NAME(dgemm)("N", "T",
                    &(C->nrow), &(C->ncol), &(B->ncol),
                    &alpha, B->elem, &(B->nrow), A->elem, &(A->nrow),
                    &beta, C->elem, &(C->nrow));

    return C;
}

/**
 * Integer square root.
 *      isqrt(n) = floor(sqrt(n))
 * which is the biggest positive integer such that its square is smaller or
 * equal to `n`.
 * 
 * @param n positive integer.
 * 
 * @returns integer square root of `n`.
 */
static int isqrt(const int n) {
    int x = 1, y;

    if (n < 0) {
        // TODO: error handling!
    } else if (n < 4) {
        return 1;
    } else if (n < 9) { // Only required for n <= 4.
        return 2;
    }

    while (1) {
        y = (x + n / x) / 2;
        if (x == y || x + 1 == y) {
            return x;
        }
        x = y;
    }
}

/**
 * Vector to symmetric matrix transform.
 * 
 * Given a vector `A` of length `N = n * (n - 1) / 2` a symmetric matrix `B`
 * of dimension `n x n` is created such that (0-indexed)
 * 
 *      b_{i j} = b_{j i} = a_k,    for    k = j * n + i - j * (j + 1) / 2
 * 
 * with 0 <= j < i < n. The main diagonal elements are set to `diag`.
 *
 * @example The vector `A = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)^T` with
 * `diag = -1` results in
 * 
 *                      ( -1  0  1  2  3  )
 *                      (  0 -1  4  5  6  )
 *                  B = (  1  4 -1  7  8  )
 *                      (  2  5  7 -1  9  )
 *                      (  3  6  8  9 -1  )
 *
 * @param A (in) source vector of length `N = n * (n - 1) / 2`.
 * @param B (in/out) destination matrix of dim. `n x n`.
 * @param diag value of all diagonal elements.
 * 
 * @returns passed matrix `B` or new created if `B` is NULL.
 */
mat* lvecToSym(const mat* A, const double diag, mat* B) {
    int i, j, k, l, n, N = A->nrow;
    int max_size, block_size;
    double *a, *b, *L, *U;

    if (!B) {
        n = (1 + isqrt(1 + 8 * N)) / 2;
        if (n * (n - 1) != 2 * N) {
            // TODO: error handling!
        }
        B = matrix(n, n);
    } else {
        n = B->nrow;
        if (n * (n - 1) != 2 * N || n != B->ncol) {
            // TODO: error handling!
        }
    }

    /* Copy consecutive values of `A` to lower triangular part of `B`. */
    a = A->elem;
    for (j = k = 0; j < n; ++j) {
        b = B->elem + j * n;
        b[j] = diag;
        for (i = j + 1; i < n; ++i, ++k) {
            b[i] = a[k];
        }
    }

    /* Mirror along the main diagonal, aka transposed copy of lower to upper. */
    b = B->elem;
    max_size = BLOCK_SIZE * (n / BLOCK_SIZE);
    /* Iterate over blocked columns */
    for (j = 0; j < n; j += BLOCK_SIZE) {
        /* determine height of the block (for reminding partial blocks) */
        block_size = j < max_size ? BLOCK_SIZE : (n % BLOCK_SIZE);
        /* Set diagonal block (left upper corner on the main diagonal) */
        L = b + (j * n + j); // diagonal block.
        /* Transpose block on main diagonal */
        for (l = 0; l < block_size; ++l) {
            for (k = l + 1; k < block_size; ++k) {
                L[k * n + l] = L[l * n + k];
            }
        }
        /* Iterate complete blocks below previous halve block */
        for (i = j + BLOCK_SIZE; i < n; i += BLOCK_SIZE) {
            /* determine height of the block (for reminding partial blocks) */
            block_size = i < max_size ? BLOCK_SIZE : (n % BLOCK_SIZE);
            /* Set lower (L, below diag) and upper (U, above diag) blocks */
            L = b + j * n + i;
            U = b + i * n + j;
            /* Transpose lower (L) to upper (U) complete block. */
            for (l = 0; l < BLOCK_SIZE; ++l) {
                for (k = 0; k < block_size; ++k) {
                    U[k * n + l] = L[l * n + k];
                }
            }
        }
    }

    return B;
}

/**
 * In place transformation into symmetric matrix as
 *      A <- diag(colSums(A + t(A))) - (A + t(A))
 *
 * @param A (in/out) square matrix.
 * @param workMem (out) double[n] vector as working memory. Will be storing the
 * diagonal elements of the result matrix. Must have at least length `n`.
 *
 * @returns laplace transformed input `A`.
 * 
 * @details The workMem parameter is required as working memory.
 */
mat* laplace(mat *A, double *workMem) {
    int i, j, k, l, n = A->nrow;
    int max_size = BLOCK_SIZE * (n / BLOCK_SIZE);
    int block_size;
    double *a = A->elem;
    double *L, *U;

    /* Check if A is square */
    if (n != A->ncol) {
        // TODO: error handling!
    }

    /* Ckeck if working memory is supplied. */
    if (!workMem) {
        workMem = (double*)R_alloc(n, sizeof(double));
    }

    /* init working memory */
    memset(workMem, 0, n * sizeof(double));

    /* Iterate over column slices */
    for (j = 0; j < n; j += BLOCK_SIZE) {
        /* determine height of the block (for reminding partial blocks) */
        block_size = j < max_size ? BLOCK_SIZE : (n % BLOCK_SIZE);
        /* Set diagonal block (left upper corner on the main diagonal) */
        L = a + (j * n + j); // D block in description.
        /* Transpose block on main diagonal */
        for (l = 0; l < block_size; ++l) {
            for (k = l + 1; k < block_size; ++k) {
                workMem[j + l] += L[l * n + k]
                                = L[k * n + l]
                                = -(L[l * n + k] + L[k * n + l]);
            }
        }
        /* Iterate complete blocks below previous halve block */
        for (i = j + BLOCK_SIZE; i < n; i += BLOCK_SIZE) {
            /* determine height of the block (for reminding partial blocks) */
            block_size = i < max_size ? BLOCK_SIZE : (n % BLOCK_SIZE);
            /* Set lower (L, below diag) and upper (U, above diag) blocks */
            L = a + j * n + i;
            U = a + i * n + j;
            /* Transpose lower (L) to upper (U) complete block. */
            for (l = 0; l < BLOCK_SIZE; ++l) {
                for (k = 0; k < block_size; ++k) {
                    workMem[j + l] += L[l * n + k]
                                    = U[k * n + l]
                                    = -(L[l * n + k] + U[k * n + l]);
                }
            }
        }
    }
    /* Sum remaining upper diagonal column sum parts and set diagonal */
    for (j = 0; j < n; ++j) {
        for (i = 0; i < j; ++i) {
            workMem[j] += a[i];
        }
        a[j] = -workMem[j];
        a += n;
    }

    return A;
}

/** Cayley transformation of matrix `B` using a Skew-Symmetric matrix `A`.
 *
 *      C = (I + A)^-1 (I - A) B
 *
 * by solving the following linear equation:
 *
 *      (I + A) C = (I - A) B ==> C = (I + A)^-1 (I - A) B
 *      \_____/     \_____/
 *        IpA   C =   ImA   B
 *                  \_______/
 *        IpA   C =     Y     ==> C = IpA^-1 Y
 *
 * @param A Skew-Symmetric matrix of dimension `(n, n)`.
 * @param B Matrix of dimensions `(n, m)` with `m <= n`.
 * @param C Matrix of dimensions `(n, m)` with `m <= n`.
 * @param workMem working memory array of length at least `2 p^2 + p`.
 * @return Transformed matrix `C`.
 * @note This opperation is equivalent to the R expression:
 *      solve(diag(1, n) + A) %*% (diag(1, n) - A) %*% B
 * or
 *      solve(diag(1, n) + A, (diag(1, n) - A) %*% B)
 */
mat* cayleyTransform(mat *A, mat *B, mat *C, double *workMem) {
    int i, info, pp = A->nrow * A->nrow;
    double zero = 0.0, one = 1.0;

    // TODO: validate dimensions!!!
    // TODO: calrify a bit ^^

    if (!C) {
        C = matrix(B->nrow, B->ncol);
    } else if (B->nrow != C->nrow || B->ncol != C->ncol) {
        // TODO: error handling!
    }

    /* Allocate row permutation array used by `dgesv` */
    int *ipiv = (int*)workMem;

    /* NOTE: workMem offset, NOT ipiv offset! There may be a bit space left out
     * but the working memory is required elsewhere anyway. It's impotant to
     * have an appropriate beginning cause if may occure the case that the
     * memory addressing faily due to size differences between int and double
     * leading to an illegal double* address. */
    double *IpA = workMem + A->nrow;
    /* Create Matrix IpA = I + A (I plus A) */
    memcpy(IpA, A->elem, A->nrow * A->ncol * sizeof(double));
    for (i = 0; i < pp; i += A->nrow + 1) {
        IpA[i] += 1.0; // +1 to diagonal elements.
    }

    /* Create Matrix ImA = I - A (I minus A) */
    double *ImA = IpA + pp;
    double *a = A->elem;
    for (i = 0; i < pp; ++i) {
        ImA[i] = -a[i];
    }
    for (i = 0; i < pp; i += A->nrow + 1) {
        ImA[i] += 1.0; // +1 to diagonal elements.
    }

    /* Y as matrix-matrix product of ImA and B:
     * Y = 1 * ImA * B + 0 * Y */
    F77_CALL(dgemm)("N", "N", &(A->nrow), &(B->ncol), &(A->nrow),
                    &one, ImA, &(A->nrow), B->elem, &(B->nrow),
                    &zero, C->elem, &(C->nrow));

    /* Solve system IpA Y = C for Y (and store result in C).
     * aka. C = IpA^-1 X */
    F77_CALL(dgesv)(&(A->nrow), &(B->ncol), IpA, &(A->nrow),
                    ipiv, C->elem, &(A->nrow), &info);

    if (info) {
        error("[ cayleyTransform ] error in dgesv - info %d", info);
    }

    return C;
}
