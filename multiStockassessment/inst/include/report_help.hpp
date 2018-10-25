
#include <string>
#include <sstream>

void moveREPORT(SEXP to, SEXP from){
  SEXP names = PROTECT(R_lsInternal(from, FALSE));
  for(int i = 0; i < Rf_length(names); ++i){
    SEXP name = PROTECT(STRING_ELT(names,i));
    // SEXP val = PROTECT(Rf_findVarInFrame(to, Rf_install(CHAR(name))));
    // if (val == R_UnboundValue) {
    SEXP valfrom = PROTECT(Rf_findVarInFrame(from, Rf_install(CHAR(name))));
    Rf_defineVar(Rf_install(CHAR(name)),valfrom,to);
    //UNPROTECT(1);
    //}
    UNPROTECT(2);
  }
  UNPROTECT(1);
  return;
}

template<class Type>
void moveADREPORT(objective_function<Type>* from, objective_function<Type>* to,int stock){  
  int nStart = 0;
  for(int i = 0; i < from->reportvector.names.size(); ++i){
    vector<int> n = from->reportvector.namedim(i);
    vector<Type> res = from->reportvector.result.segment(nStart,n.prod());
    nStart += n.prod();
    std::string nam("");
    if(stock >= 0){
      nam.append("SAM_");
      std::ostringstream s;
      s << stock;
      nam.append(s.str());
      nam.append("_");
    }
    nam.append(from->reportvector.names(i));
    to->reportvector.push(res,strdup(nam.data()));
  }
  return;
}


template<class Type>
struct oftmp :
  public objective_function<Type> {

  oftmp() :
    objective_function<Type>(PROTECT(Rf_allocVector(VECSXP,0)),
			     PROTECT(Rf_allocVector(VECSXP,0)),
			     PROTECT(allocSExp(ENVSXP)))
  {};

  oftmp(bool do_simulate):
    objective_function<Type>(PROTECT(Rf_allocVector(VECSXP,0)),
			     PROTECT(Rf_allocVector(VECSXP,0)),
			     PROTECT(allocSExp(ENVSXP)))
  {
    this->set_simulate(do_simulate);
  };
  
  Type operator()(){
    return Type();
  }
  
  ~oftmp(){
    UNPROTECT(3);
  }
};

template<class Type>
struct ofall :
  public objective_function<Type> {

  int nStocks;
  
  ofall(int nStocks_) :
    objective_function<Type>(PROTECT(Rf_allocVector(VECSXP,0)),
			     PROTECT(Rf_allocVector(VECSXP,0)),
			     PROTECT(allocSExp(ENVSXP))),
    nStocks(nStocks_)
  {};

  ~ofall(){
    UNPROTECT(3);
  }

  Type operator()(){
    return Type();
  }

  void addToReport(SEXP rep, int stock){
    if(!isDouble<Type>::value) return;
    SEXP names = PROTECT(R_lsInternal(rep, FALSE));
    if(names == NILSXP) return;
    for(int i = 0; i < Rf_length(names); ++i){
      SEXP name = PROTECT(STRING_ELT(names,i));
      SEXP val = PROTECT(Rf_findVarInFrame(this->report, Rf_install(CHAR(name))));
      if (val == R_UnboundValue) {
	SEXP vec = PROTECT(Rf_allocVector(VECSXP,nStocks));
	SEXP newval = PROTECT(Rf_findVarInFrame(rep,Rf_install(CHAR(name))));
	SET_VECTOR_ELT(vec,stock,newval);
	Rf_defineVar(Rf_install(CHAR(name)),vec,this->report);
	UNPROTECT(2);
      }else{
	SEXP newval = PROTECT(Rf_findVarInFrame(rep,Rf_install(CHAR(name))));
	SET_VECTOR_ELT(val,stock,newval);
	UNPROTECT(1);
      }
      UNPROTECT(2);
    }
    UNPROTECT(1);
    return;
  }
  
};
