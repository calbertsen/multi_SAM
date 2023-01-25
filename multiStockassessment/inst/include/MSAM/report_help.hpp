
#include <string>
#include <sstream>

void moveREPORT(SEXP to, SEXP from)SOURCE({
   SEXP names = PROTECT(R_lsInternal(from, FALSE));
    if(names == NILSXP){
      return;
    }
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
				PROTECT(allocSExp(ENVSXP)))
       {};
       );

SOURCE(
       template<class Type>
       oftmp<Type>::oftmp(bool do_simulate):
       objective_function<Type>(PROTECT(Rf_allocVector(VECSXP,0)),
				PROTECT(Rf_allocVector(VECSXP,0)),
				PROTECT(allocSExp(ENVSXP)))
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
				PROTECT(allocSExp(ENVSXP))),
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
	 SEXP names = PROTECT(R_lsInternal(rep, FALSE));
	 if(names == NILSXP){
	   return;
	 }
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
       );

MSM_SPECIALIZATION(struct ofall<double>);
MSM_SPECIALIZATION(struct ofall<TMBad::ad_aug>);
