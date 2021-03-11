#include "cve.h"

/**
 * Calls a R function passed to the algorithm and supplied intermediate
 * optimization values for logging the optimization progress.
 * The supplied parameters to the logger functions are as follows:
 * - attempt: Attempts counter.
 * - iter: Current iter staring with 0 as initial iter.
 * - L: Per X_i to X_j pair loss.
 * - V: Current estimated SDR null space basis.
 * - tau: Step-size.
 * - err: Error \eqn{|| V V^T - V_{tau} V_{tau}^T ||}.
 *
 * @param logger Pointer to a SEXP R object representing an R function.
 * @param env Pointer to a SEXP R object representing an R environment.
 * @param attempt counter of attempts.
 * @param iter optimization iteration counter.
 * @param L Pointer to a SEXP R object representing an R environment.
 * @param V Pointer memory area of size `nrowV * ncolV` storing `V`.
 * @param G Pointer memory area of size `nrowG * ncolG` storing `G`.
 * @param loss Current loss L(V).
 * @param err Error for break condition (0.0 before first iteration).
 * @param tau Current step-size.
 */
void callLogger(SEXP logger, SEXP env,
                const int attempt, const int iter,
                const mat* L, const mat* V, const mat* G,
                const double loss, const double err, const double tau) {
    /* Create R objects to be passed to R logger function. */
    // Attempt is converted from 0-indexed to 1-indexed as R index.
    SEXP r_attempt = PROTECT(ScalarInteger(attempt + 1));
    SEXP r_iter = PROTECT(ScalarInteger(iter + 1));

    /* Create R representations of L, V and G */
    SEXP r_L = PROTECT(allocMatrix(REALSXP, L->nrow, L->ncol));
    SEXP r_V = PROTECT(allocMatrix(REALSXP, V->nrow, V->ncol));
    SEXP r_G = PROTECT(allocMatrix(REALSXP, G->nrow, G->ncol));
    /* Copy data to R objects */
    memcpy(REAL(r_L), L->elem, L->nrow * L->ncol * sizeof(double));
    memcpy(REAL(r_V), V->elem, V->nrow * V->ncol * sizeof(double));
    memcpy(REAL(r_G), G->elem, G->nrow * G->ncol * sizeof(double));

    /* Build data list passed to logger */
    SEXP data = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(data, 0, r_L);
    SET_VECTOR_ELT(data, 1, r_V);
    SET_VECTOR_ELT(data, 2, r_G);
    SET_VECTOR_ELT(data, 3, PROTECT(ScalarReal(loss)));
    SET_VECTOR_ELT(data, 4, PROTECT(ScalarReal(err < 0.0 ? NA_REAL : err)));
    SET_VECTOR_ELT(data, 5, PROTECT(ScalarReal(tau)));
    SEXP names = PROTECT(allocVector(STRSXP, 6));
    SET_STRING_ELT(names, 0, mkChar("L"));
    SET_STRING_ELT(names, 1, mkChar("V"));
    SET_STRING_ELT(names, 2, mkChar("G"));
    SET_STRING_ELT(names, 3, mkChar("loss"));
    SET_STRING_ELT(names, 4, mkChar("err"));
    SET_STRING_ELT(names, 5, mkChar("tau"));
    setAttrib(data, R_NamesSymbol, names);

    /* Create logger function call as R language expression. */
    SEXP loggerCall = PROTECT(lang4(logger, r_attempt, r_iter, data));

    /* Evaluate the logger function call expression. */
    eval(loggerCall, env);

    /* Unlock created R objects. */
    UNPROTECT(11);
}
