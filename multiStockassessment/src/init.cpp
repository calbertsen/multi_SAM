#define WITH_LIBTMB
#include <TMB.hpp>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


extern "C" {

#define CALLDEF(name,n) {#name, (DL_FUNC) &name, n}
  
  static const
  R_CallMethodDef callMethods[] = {

    // TMB
    #ifdef TMB_CALLDEFS
    TMB_CALLDEFS,
    #else
    CALLDEF(MakeADFunObject, 4),
    CALLDEF(InfoADFunObject, 1),
    CALLDEF(EvalADFunObject, 3),
    CALLDEF(MakeDoubleFunObject, 3),
    CALLDEF(EvalDoubleFunObject, 3),
    CALLDEF(getParameterOrder, 3),
    CALLDEF(MakeADGradObject, 3),
    CALLDEF(MakeADHessObject2, 4),
    CALLDEF(usingAtomics, 0),
    CALLDEF(TMBconfig, 2),
    #endif
    
    // CALLDEF(rmvnorm,3),
    // CALLDEF(rmvlnorm,3),
    // CALLDEF(rlogistic,3),
    // CALLDEF(rcategory,2),
    // CALLDEF(rmvgaussmix,4),
    // CALLDEF(rmvt,4),
    // CALLDEF(rgig,4),
    // CALLDEF(ll2n,2),
    // CALLDEF(n2ll,2),
    // CALLDEF(idtcrwVarMat,5),
    // CALLDEF(stepLength,5),
    // CALLDEF(bearing,5),

    {NULL,NULL,0}
  };

  void R_init_multiStockassessment(DllInfo *info)
  {
    /* Register the .C and .Call routines.
       No .Fortran() or .External() routines,
       so pass those arrays as NULL.
    */
    R_registerRoutines(info,
		       NULL, callMethods,
		       NULL, NULL);
    R_useDynamicSymbols(info, (Rboolean)FALSE);
  }


}
