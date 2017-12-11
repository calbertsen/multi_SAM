
#include <string>
#include <sstream>

void addMissingVars(SEXP to, SEXP from){
  SEXP names = R_lsInternal(from, FALSE);
  for(int i = 0; i < Rf_length(names); ++i){
    SEXP name = STRING_ELT(names,i);
    SEXP val = Rf_findVarInFrame(to, Rf_install(CHAR(name)));
    if (val == R_UnboundValue) {
	SEXP valfrom = Rf_findVarInFrame(from, Rf_install(CHAR(name)));
	Rf_defineVar(Rf_install(CHAR(name)),valfrom,to);
    }
  }
  return;
}

template<class Type>
void moveADREPORT(objective_function<Type>* from, objective_function<Type>* to,int stock){  
  int nStart = 0;
  for(int i = 0; i < from->reportvector.names.size(); ++i){
    int n = from->reportvector.namelength(i);
    vector<Type> res = from->reportvector.result.segment(nStart,n);
    nStart += n;
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
    objective_function<Type>(Rf_allocVector(VECSXP,0),
			     Rf_allocVector(VECSXP,0),
			     allocSExp(ENVSXP))
  {};
  Type operator()(){
    return Type();
  }
};

template<class Type>
struct ofall :
  public objective_function<Type> {

  int nStocks;
  
  ofall(int nStocks_) :
    objective_function<Type>(Rf_allocVector(VECSXP,0),
			     Rf_allocVector(VECSXP,0),
			     allocSExp(ENVSXP)),
    nStocks(nStocks_)
  {};
  Type operator()(){
    return Type();
  }

  void addToReport(SEXP rep, int stock){
    if(!isDouble<Type>::value) return;
    SEXP names = R_lsInternal(rep, FALSE);
    if(names == NILSXP) return;
    for(int i = 0; i < Rf_length(names); ++i){
      SEXP name = STRING_ELT(names,i);
      SEXP val = Rf_findVarInFrame(this->report, Rf_install(CHAR(name)));
      if (val == R_UnboundValue) {
	SEXP vec = Rf_allocVector(VECSXP,nStocks);
	SEXP newval = Rf_findVarInFrame(rep,Rf_install(CHAR(name)));
	SET_VECTOR_ELT(vec,stock,newval);
	Rf_defineVar(Rf_install(CHAR(name)),vec,this->report);
      }else{
	SEXP newval = Rf_findVarInFrame(rep,Rf_install(CHAR(name)));
	SET_VECTOR_ELT(val,stock,newval);
      }
    }
    return;
  }
  
};
