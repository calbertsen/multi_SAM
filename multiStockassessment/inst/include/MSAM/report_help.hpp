
#include <string>
#include <sstream>

SEXP MSAM_R_ls(SEXP env, Rboolean all) SOURCE({
  SEXP c = PROTECT(lang4(Rf_install("ls"),
			 env,
			 Rf_ScalarLogical(all),
			 Rf_ScalarLogical(TRUE)
			 ));
  
  SET_TAG(CDR(c),    Rf_install("envir"));
  SET_TAG(CDDR(c),   Rf_install("all.names"));
  SET_TAG(CDDDR(c),   Rf_install("sorted"));

  SEXP r = PROTECT(Rf_eval(c, R_GlobalEnv));
  UNPROTECT(2);
  return r;
  });

#if R_VERSION >= R_Version(4,5,0)
SEXP MSAM_R_findVar(SEXP env, SEXP name) SOURCE({    
    // No need to test with R_existsVarInFrame, R_getVar will give an error if not found
    // WRONG: we do need to test because the code after expects R_UnboundValue if it wasn't found
    SEXP sym;
    if (TYPEOF(name) == SYMSXP) {
      sym = name;
    } else {
      sym = Rf_install(CHAR(Rf_asChar(name)));
    }
    if(!R_existsVarInFrame(env,sym))
      return R_UnboundValue;
    return R_getVar(sym, env, FALSE);
  });
#else
SEXP MSAM_R_findVar(SEXP env, SEXP name) SOURCE({
    SEXP sym;
    if (TYPEOF(name) == SYMSXP) {
      sym = name;
    } else {
      sym = Rf_install(CHAR(Rf_asChar(name)));
    }
    return Rf_findVarInFrame(env, sym);
  });
#endif

#if R_VERSION >= R_Version(4,1,0)
SEXP MSAM_R_NewEnv() SOURCE({
    return R_NewEnv(R_GlobalEnv, FALSE, 0);
  });
#else
SEXP MSAM_R_NewEnv() SOURCE({
    return eval(lang1(Rf_install("new.env")), R_GlobalEnv);
  });
#endif

void moveREPORT(SEXP to, SEXP from)SOURCE({
    SEXP names = PROTECT(MSAM_R_ls(from, FALSE));
    if(names == NILSXP){
      return;
    }
  for(int i = 0; i < Rf_length(names); ++i){
    SEXP name = PROTECT(STRING_ELT(names,i));
    // SEXP val = PROTECT(Rf_findVarInFrame(to, Rf_install(CHAR(name))));
    // if (val == R_UnboundValue) {
    SEXP sym;
    if (TYPEOF(name) == SYMSXP) {
      sym = name;
    } else {
      sym = Rf_install(CHAR(Rf_asChar(name)));
    }
    SEXP valfrom = PROTECT(MSAM_R_findVar(from, sym));
    Rf_defineVar(sym,valfrom,to);
    //UNPROTECT(1);
    //}
    UNPROTECT(2);
  }
   UNPROTECT(1);
  return;
})

  
template<class Type>
std::vector<Type> segment(std::vector<Type>& x, int i, int n)SOURCE({
  return std::vector<Type>(x.begin() + i, x.begin() + i + n);
  })

  MSM_SPECIALIZATION(std::vector<double> segment(std::vector<double>&, int, int));
MSM_SPECIALIZATION(std::vector<TMBad::ad_aug> segment(std::vector<TMBad::ad_aug>&, int, int));
			  

template<class Type>
void moveADREPORT(objective_function<Type>* from, objective_function<Type>* to,int stock)SOURCE({  
  int nStart = 0;
  for(int i = 0; i < (int)from->reportvector.names.size(); ++i){
    vector<int> n = from->reportvector.namedim[i];
    std::vector<Type> rtmp = segment(from->reportvector.result,nStart,n.prod());
    vector<Type> res(rtmp);
    nStart += n.prod();
    std::string nam("");
    if(stock >= 0){
      nam.append("SAM_");
      std::ostringstream s;
      s << stock;
      nam.append(s.str());
      nam.append("_");
    }
    nam.append(from->reportvector.names[i]);
    to->reportvector.push(res,strdup(nam.data()));
  }
  return;
  })

MSM_SPECIALIZATION(void moveADREPORT(objective_function<double>* from, objective_function<double>* to,int stock));
MSM_SPECIALIZATION(void moveADREPORT(objective_function<TMBad::ad_aug>* from, objective_function<TMBad::ad_aug>* to,int stock));

HEADER(
template<class Type>
struct oftmp :
  public objective_function<Type> {

  oftmp();

  oftmp(bool do_simulate);
  
  Type operator()();
  
  ~oftmp();
};
       )

SOURCE(
       template<class Type>
       oftmp<Type>::oftmp() :
       objective_function<Type>(PROTECT(Rf_allocVector(VECSXP,0)),
				PROTECT(Rf_allocVector(VECSXP,0)),
				PROTECT(MSAM_R_NewEnv()))
       {};
       );

SOURCE(
       template<class Type>
       oftmp<Type>::oftmp(bool do_simulate):
       objective_function<Type>(PROTECT(Rf_allocVector(VECSXP,0)),
				PROTECT(Rf_allocVector(VECSXP,0)),
				PROTECT(MSAM_R_NewEnv()))
       {
	 this->set_simulate(do_simulate);
       };
       );

SOURCE(
       template<class Type>
       Type oftmp<Type>::operator()(){
	 return Type();
       });
  
SOURCE(
       template<class Type>
       oftmp<Type>::~oftmp(){
	 UNPROTECT(3);
       });


MSM_SPECIALIZATION(struct oftmp<double>);
MSM_SPECIALIZATION(struct oftmp<TMBad::ad_aug>);


HEADER(
       template<class Type>
       struct ofall :
       public objective_function<Type> {

	 int nStocks;
  
	 ofall(int nStocks_);
	 ~ofall();

	 Type operator()();

	 void addToReport(SEXP rep, int stock);  
       };
       );

SOURCE(
       template<class Type>
       ofall<Type>::ofall(int nStocks_) :
       objective_function<Type>(PROTECT(Rf_allocVector(VECSXP,0)),
				PROTECT(Rf_allocVector(VECSXP,0)),
				PROTECT(MSAM_R_NewEnv())),
       nStocks(nStocks_)
       {};
       );

SOURCE(
       template<class Type>
       ofall<Type>::~ofall(){
	 UNPROTECT(3);
       }
       );

SOURCE(
       template<class Type>
       Type ofall<Type>::operator()(){
	 return Type();
       }
       );

SOURCE(
       template<class Type>
       void ofall<Type>::addToReport(SEXP rep, int stock){
	 if(!isDouble<Type>::value) return;
	 SEXP names = PROTECT(MSAM_R_ls(rep, FALSE));
	 if(names == NILSXP){
	   return;
	 }
	 for(int i = 0; i < Rf_length(names); ++i){
	   SEXP name = PROTECT(STRING_ELT(names,i));
	   SEXP sym;
	   if (TYPEOF(name) == SYMSXP) {
	     sym = PROTECT(name);
	   } else {
	     sym = PROTECT(Rf_install(CHAR(Rf_asChar(name))));
	   }

	   SEXP val = PROTECT(MSAM_R_findVar(this->report, sym));
	   if (val == R_UnboundValue) {
	     SEXP vec = PROTECT(Rf_allocVector(VECSXP,nStocks));
	     SEXP newval = PROTECT(MSAM_R_findVar(rep,sym));
	     SET_VECTOR_ELT(vec,stock,newval);
	     Rf_defineVar(sym,vec,this->report);
	     UNPROTECT(2);
	   }else{
	     SEXP newval = PROTECT(MSAM_R_findVar(rep,sym));
	     SET_VECTOR_ELT(val,stock,newval);
	     UNPROTECT(1);
	   }
	   UNPROTECT(3);
	 }
	 UNPROTECT(1);
	 return;
       }
       );

MSM_SPECIALIZATION(struct ofall<double>);
MSM_SPECIALIZATION(struct ofall<TMBad::ad_aug>);
