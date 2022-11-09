#ifndef MSAM_DEPENDS
#define MSAM_DEPENDS(x) ;
#endif

#ifndef SAM_DEPENDS
#define SAM_DEPENDS(x) ;
#endif


#ifdef WITH_MSAM_LIB		/* Only get headers */
#undef SOURCE
#define SOURCE(...) ;
#else
#undef SOURCE
#define SOURCE(...) __VA_ARGS__;
#endif

#ifdef MSAM_COMPILEUNITS
#ifdef WITH_MSAM_LIB
#undef DEFARG
#define DEFARG(...) __VA_ARGS__
#undef MSAM_SPECIALIZATION
#define MSAM_SPECIALIZATION(...) extern template __VA_ARGS__
#undef HEADER
#define HEADER(...) __VA_ARGS__;
#else 
#undef DEFARG
#define DEFARG(...)
#undef MSAM_SPECIALIZATION
#define MSAM_SPECIALIZATION(...) template __VA_ARGS__
#undef HEADER
#define HEADER(...);
#endif
#else  /* Use as header only library */
#undef DEFARG
#define DEFARG(...) __VA_ARGS__
#undef MSAM_SPECIALIZATION
#define MSAM_SPECIALIZATION(...) ;
#undef HEADER
#define HEADER(...) __VA_ARGS__;
#endif



#ifndef MSAM_ASSERT
#define MSAM_ASSERT(x,msg) if(!(x)) Rf_error(msg);
#endif
