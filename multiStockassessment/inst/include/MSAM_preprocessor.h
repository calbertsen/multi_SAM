#ifndef MSM_DEPENDS
#define MSM_DEPENDS(x) ;
#endif

#ifndef SAM_DEPENDS
#define SAM_DEPENDS(x) ;
#endif


#ifdef WITH_MSM_LIB		/* Only get headers */
#undef SOURCE
#define SOURCE(...) ;
#else
#undef SOURCE
#define SOURCE(...) __VA_ARGS__;
#endif

#ifdef MSM_COMPILEUNITS
#ifdef WITH_MSM_LIB
#undef DEFARG
#define DEFARG(...) __VA_ARGS__
#undef MSM_SPECIALIZATION
#define MSM_SPECIALIZATION(...) extern template __VA_ARGS__
#undef HEADER
#define HEADER(...) __VA_ARGS__;
#else 
#undef DEFARG
#define DEFARG(...)
#undef MSM_SPECIALIZATION
#define MSM_SPECIALIZATION(...) template __VA_ARGS__
#undef HEADER
#define HEADER(...);
#endif
#else  /* Use as header only library */
#undef DEFARG
#define DEFARG(...) __VA_ARGS__
#undef MSM_SPECIALIZATION
#define MSM_SPECIALIZATION(...) ;
#undef HEADER
#define HEADER(...) __VA_ARGS__;
#endif



#ifndef MSM_ASSERT
#define MSM_ASSERT(x,msg) if(!(x)) Rf_error(msg);
#endif
