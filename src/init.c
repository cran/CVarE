#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP cve_export(SEXP X, SEXP Fy, SEXP k, SEXP h,
                       SEXP method,
                       SEXP V, // initial
                       SEXP momentum, SEXP tau, SEXP tol,
                       SEXP slack, SEXP gamma,
                       SEXP maxIter, SEXP attempts,
                       SEXP logger, SEXP loggerEnv);

static const R_CallMethodDef CallEntries[] = {
    {"cve_export", (DL_FUNC) &cve_export, 15},
    {NULL, NULL, 0}
};

/* Restrict C entry points to registered routines. */
void R_init_CVE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
