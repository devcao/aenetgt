#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP aenetgt_CovYiYjgibbs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP aenetgt_EYgibbs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP aenetgt_EYiYjgibbs_slow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP aenetgt_llj_array(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP aenetgt_logistic_enet(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"aenetgt_CovYiYjgibbs",    (DL_FUNC) &aenetgt_CovYiYjgibbs,    10},
    {"aenetgt_EYgibbs",         (DL_FUNC) &aenetgt_EYgibbs,          8},
    {"aenetgt_EYiYjgibbs_slow", (DL_FUNC) &aenetgt_EYiYjgibbs_slow,  8},
    {"aenetgt_llj_array",       (DL_FUNC) &aenetgt_llj_array,        8},
    {"aenetgt_logistic_enet",   (DL_FUNC) &aenetgt_logistic_enet,    7},
    {NULL, NULL, 0}
};

void R_init_aenetgt(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

