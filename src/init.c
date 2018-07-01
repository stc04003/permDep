// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdlib.h> // for NULL

/* .C calls */
extern void kendallTrun(void *, void *, void *, void *, void *);
extern void kendallTrunWgt(void *, void *, void *, void *, void *, void *, void *);
extern void mysampleC(void *, void *, void *, void *);

/* .Call calls */
extern SEXP agfit4(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"kendallTrun",    (DL_FUNC) &kendallTrun,    5},
    {"kendallTrunWgt", (DL_FUNC) &kendallTrunWgt, 7},
    {"mysampleC",      (DL_FUNC) &mysampleC,      4},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"agfit4", (DL_FUNC) &agfit4, 13},
    {NULL, NULL, 0}
};

void R_init_permDep(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
