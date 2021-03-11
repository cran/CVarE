/* Include Guard */
#ifndef CVE_INCLUDE_GUARD_H_
#define CVE_INCLUDE_GUARD_H_

#include <string.h> // `mem*` functions.
#include <math.h>   // sqrt, exp, ...

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#define CVE_MEM_CHUNK_SIZE 2032
#define CVE_MEM_CHUNK_SMALL 1016
#define BLOCK_SIZE 8

/**
 * @struct Matrix of dimensions `nrow x ncol`.
 */
typedef struct matrix {
    int nrow;     /**< Number of rows */
    int ncol;     /**< Number of columns */
    void *origin; /**< Reference to origin, see `asMat()`. */
    double *elem; /**< Column-major array of matrix elements. */
} mat;

typedef enum {
    simple,
    weighted
} method;

typedef enum {
    gauss
} kernel;


void cve(const mat *X, const mat *Y, const double h,
         const unsigned int method,
         const double momentum,
         const double tau_init, const double tol_init,
         const double slack, const double gamma,
         const int maxIter, const int attempts,
         mat *V, mat *L,
         SEXP logger, SEXP loggerEnv);

void callLogger(SEXP logger, SEXP env,
                const int attempt, const int iter,
                const mat* L, const mat* V, const mat* G,
                const double loss, const double err, const double tau);

/******************************************************************************/
/**                                rStiefel.c                                **/
/******************************************************************************/
/* Random element from Stiefel manifold. */
mat* rStiefel(const int p, const int q, mat *V, double *workMem);

/******************************************************************************/
/**                                 matrix.c                                 **/
/******************************************************************************/
/* Create and Copy matrices */
mat* matrix(const int nrow, const int ncol);
mat* zero(const int nrow, const int ncol);
mat* copy(mat *src, mat *dest);
/* Matrix to scalar */
double sum(const mat *A);
double mean(const mat *A);
double squareSum(const mat* A);
double dot(const mat *A, const char op, const mat *B);
double dist(const mat *A, const mat *B);
/* Matrix to vector (`ncol == 1` matrices, aka row vectors) */
mat* rowSums(const mat *A, mat *sums);
mat* colSums(const mat *A, mat *sums);
mat* rowDiffSquareSums(const mat *X, mat *lvecD);
/* Matrix and scalar to Matrix */
mat* elemApply(mat *A, const char op, const double scalar, mat *B);
/* Matrix and vector to Matrix */
mat* colApply(mat *A, const char op, mat *B, mat *C);
/* Matrix and Matrix to Matrix */
mat* lincomb(const double alpha, const mat *A, const double beta, mat *B);
mat* matrixprod(const double alpha, const mat *A, const mat *B,
                double beta, mat* C);
mat* crossprod(const double alpha, const mat *A, const mat *B,
               double beta, mat* C);
mat* hadamard(const double alpha, const mat* A, const mat *B,
              double beta, mat* C);
mat* skew(const double alpha, mat* A, mat *B, double beta, mat* C);
/* Matrix Transformations */
mat* cayleyTransform(mat *A, mat *B, mat *C, double *workMem);
mat* laplace(mat *A, double *workMem);
mat* lvecToSym(const mat* A, const double diag, mat* B);

/******************************************************************************/
/**                            cve_subroutines.c                             **/
/******************************************************************************/
/* CVE specific sub-routines */
mat* adjacence(const mat *vec_L, const mat *vec_Y, const mat *vec_y1,
               const mat *mat_D, const mat *mat_W, kernel kernel,
               mat *mat_S);
mat* applyKernel(const mat* A, const double h, kernel kernel, mat* B);

#endif /* CVE_INCLUDE_GUARD_H_ */
