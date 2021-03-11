#include "cve.h"

/**
 * Converts a `SEXP` (S EXPression) into a matrix struct `mat`.
 * 
 * @param S source struct to be converted.
 *
 * @details Reuses the memory area of the SEXP object, therefore manipulation
 *      of the returned matrix works in place of the SEXP object. In addition,
 *      a reference to the original SEXP is stored and will be retrieved from
 *      `asSEXP()` if the matrix was created through this function.
 */
static mat* asMat(SEXP S) {
    // TODO: error checking and conversion
    mat* M = (mat*)R_alloc(1, sizeof(mat));
    if (isMatrix(S)) {
        M->nrow = (int)nrows(S);
        M->ncol = (int)ncols(S);
    } else {
        M->nrow = (int)length(S);
        M->ncol = 1;
    }
    M->origin = S;
    M->elem = REAL(S);
    return M;
}

SEXP cve_export(SEXP X, SEXP Fy, SEXP k, SEXP h,
                SEXP method,
                SEXP V, // initial
                SEXP momentum, SEXP tau, SEXP tol,
                SEXP slack, SEXP gamma,
                SEXP maxIter, SEXP attempts,
                SEXP logger, SEXP loggerEnv) {
    /* Handle logger parameter, set to NULL pointer if not a function. */
    if (!(isFunction(logger) && isEnvironment(loggerEnv))) {
        logger = (void*)0;
    }

    /* Get dimensions. */
    int n = nrows(X);
    int p = ncols(X);
    int q = p - asInteger(k);

    /* Convert types if needed. */
    // TODO: implement! (or leave in calling R code?)

    /* Create output list. */
    SEXP Vout = PROTECT(allocMatrix(REALSXP, p, q));
    SEXP Lout = PROTECT(allocMatrix(REALSXP, n, ncols(Fy)));

    /* Check `attempts`, if not positive use passed values of `V` as
     * optimization start value without further attempts.
     * Therefor, copy from `V` to `Vout`. */
    if (asInteger(attempts) < 1) {
        // TODO: Check for 
        memcpy(REAL(Vout), REAL(V), p * q * sizeof(double));
    }

    /* Call CVE */
    cve(asMat(X), asMat(Fy), asReal(h),
        asInteger(method),
        asReal(momentum), asReal(tau), asReal(tol),
        asReal(slack), asReal(gamma),
        asInteger(maxIter), asInteger(attempts),
        asMat(Vout), asMat(Lout),
        logger, loggerEnv);

    /* Build output list object with names "V", "L" */
    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, Vout);
    SET_VECTOR_ELT(out, 1, Lout);
    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("V"));
    SET_STRING_ELT(names, 1, mkChar("L"));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(4);
    return out;
}
