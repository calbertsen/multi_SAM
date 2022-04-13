
// template<class Type>
// struct TOTAL_Z_y {

//   matrix<Type> F;		// One col per stock
//   matrix<Type> M;
//   matrix<Type> N;

//   template<class T>
//   TOTAL_Z_y(const matrix<T>& F_,
// 	    const matrix<T>& M_,
// 	    const matrix<T>& N_) :
//     F(F_), M(M_), N(N_) {}

//   Type operator()(const vector<Type> &x) {
//     // Length of x should be 2 * length of F     
//     int nage = x.size()/2;
//     int nstock = F.cols();
//     Type nll = 0.0;
//     for(int a = 0; a < nage; ++a){
//       Type vN = 0.0;
//       Type vC = 0.0;
//       Type Nat = 0.0;
//       for(int s = 0; s < nstock; ++s){
// 	Type Za = F(a,s) + M(a,s);
// 	vN += exp(-Za) * N(a,s);
// 	vC += F(a,s) / Za * (1.0 - exp(-Za)) * N(a,s);
// 	Nat += N(a,s);
//       }
//       Type Fat = x(a);
//       Type Zat = x(a + nage);
//       Type tmp1 = vN - exp(-Zat) * Nat;
//       Type tmp2 = vC - Fat / Zat * (1.0 - exp(-Zat)) * Nat;       
//       nll += tmp1 * tmp1 + tmp2 * tmp2;
//     }
//     return nll;
//   }

//   template<class T>
//   T toLogFbar(const vector<T> &x){
//     int nage = x.size()/2;
//     T logFsum = R_NegInf;
//     for(int a = 0; a < nage; ++a)
//       logFsum = logspace_add2(logFsum, x(a));
//     return logFsum - log(T(nage));
//   }
// };


// template<class Type>
// vector<Type> totalLogFbar(vector<confSet> confA,
// 			  vector<dataSet<Type> > datA,		       		       
// 		       cmoe_matrix<Type> logF,
// 		       cmoe_matrix<Type> logN,
// 		       int minYearAll,
// 		       int maxYearAll
// 		       ){

//   int nYear = maxYearAll - minYearAll + 1;
//   int A1 = confA(0).fbarRange(0);
//   int A2 = confA(0).fbarRange(1);
//   int nAge = A2 - A1 + 1;
//   int nStock = confA.size();

//   vector<Type> logFbar(nYear);
//   logFbar.setConstant(R_NaReal);
//   // Loop through years
//   for(int yall = 0; yall < maxYearAll - minYearAll + 1; ++yall){
//     // Declare matrices
//     matrix<Type> F(nAge, nStock); F.setZero();
//     matrix<Type> M(nAge, nStock); M.setZero();
//     matrix<Type> N(nAge, nStock); N.setConstant(R_NaReal);
//     // Fill matrices
//     bool skip = false;
//     for(int s = 0; s < nStock; ++s){
//       matrix<Type> logFa = logF.col(s);
//       matrix<Type> logNa = logN.col(s);
//       matrix<Type> natMor = datA(s).natMor;
//       int y = yall - CppAD::Integer(datA(s).years(0) - minYearAll);      
//       if(y >= 0 && y < datA(s).noYears){	
// 	for(int aall = A1; aall <= A2; ++aall){
// 	  int a = aall - confA(s).minAge;
// 	  if(a > 0 && a < logNa.rows()){
// 	    if(confA(s).keyLogFsta(0,a) > 0)
// 	      F(a,s) = exp(logFa(confA(s).keyLogFsta(0,a),y));
// 	    M(a,s) = natMor(y,a);
// 	    N(a,s) = exp(logNa(a,y));
// 	  }else if(a >= logNa.rows()){
// 	    int au = logNa.rows() - 1;
// 	    if(confA(s).keyLogFsta(0,au) > 0)
// 	      F(a,s) = exp(logFa(confA(s).keyLogFsta(0,au),y));
// 	    M(a,s) = natMor(y,au);
// 	    N(a,s) = exp(logNa(au,y));
// 	  }
// 	}	
//       }else{
// 	skip = true;
//       }
//     }
//     if(!skip){
//       TOTAL_Z_y<TMBad::ad_aug> Fun(F,M,N);
//       vector<Type> xx(nAge * 2);
//       xx.setZero();
//       xx = newton::Newton(Fun,xx);
//       logFbar(yall) = Fun.toLogFbar(xx);
//     }
//   }
//   return logFbar;
// }


template<class Type>
void getTotals(vector<dataSet<Type> > datA,
	       vector<confSet> confA,
	       vector<paraSet<Type> > parA,
	       cmoe_matrix<Type> logF,
	       cmoe_matrix<Type> logN,
	       int minYearAll,
	       int maxYearAll,
	       int minAgeAll,
	       int maxAgeAll,
	       objective_function<Type> *of){


  // Calculate values to report
  int nYear = maxYearAll - minYearAll + 1;
  int nAge = maxAgeAll - minAgeAll + 1;
  // int nStock = datA.size();
  vector<Type> total_logssb(nYear);
  total_logssb.setConstant(R_NegInf);
  vector<Type> total_logfsb(nYear);
  total_logfsb.setConstant(R_NegInf);
  vector<Type> total_logCatch(nYear);
  total_logCatch.setConstant(R_NegInf);
  matrix<Type> total_logCatAge(nAge,nYear);
  total_logCatAge.setConstant(R_NegInf);
  vector<Type> total_logLand(nYear);
  total_logLand.setConstant(R_NegInf);
  vector<Type> total_logtsb(nYear);
  total_logtsb.setConstant(R_NegInf);
  vector<Type> total_logR(nYear);  
  total_logR.setConstant(R_NegInf);
  vector<Type> total_logfbar(nYear);
  total_logfbar.setConstant(R_NegInf);
  // vector<Type> total_logfbar = totalLogFbar(confA, datA, logF, logN, minYearAll, maxYearAll);
  vector<Type> total_logfbarL(nYear);
  total_logfbarL.setConstant(R_NegInf);

  vector<Type> total_nStocks(nYear);
  total_nStocks.setConstant(R_NegInf);

  // matrix<Type> totalNay(nAge,nYear);
  // vector<matrix<Type> > totalMay(datA.size());
  // vector<matrix<Type> > totalFay(datA.size());
  // for(int s = 0; s < datA.size(); ++s){
    
  // }
    
  for(int s = 0; s < datA.size(); ++s){
    dataSet<Type> dat = datA(s);
    confSet conf = confA(s);
    paraSet<Type> par = parA(s);

    array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
    logNa = logN.col(s);
    array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
    logFa = logF.col(s);
    
    // Calculate values to report
    
    vector<Type> S_logssb = ssbFun(dat, conf, logNa, logFa, true);
    vector<Type> S_logfsb = log(fsbFun(dat, conf, logNa, logFa));
    vector<Type> S_logCatch = catchFun(dat, conf, logNa, logFa, true);
    matrix<Type> S_logCatAge = catchFunAge(dat, conf, logNa, logFa, true);
    vector<Type> S_logLand = log(landFun(dat, conf, logNa, logFa));
    // vector<Type> varLogCatch = varLogCatchFun(dat, conf, logN, logF, par);
    // vector<Type> varLogLand = varLogLandFun(dat, conf, logN, logF, par);
    vector<Type> S_logtsb = log(tsbFun(dat, conf, logNa));
    vector<Type> S_logR = log(rFun(logNa));
    vector<Type> S_logfbar = fbarFun(conf, logFa, true);
    vector<Type> S_logfbarL = log(landFbarFun(dat, conf, logFa));
    
    for(int yall = 0; yall < maxYearAll - minYearAll + 1; ++yall){    
      int y = yall - CppAD::Integer(dat.years(0) - minYearAll);
      if(y >= 0 && y < dat.noYears){	
	if(y < S_logssb.size())
	  total_logssb(yall) = logspace_add2(total_logssb(yall),S_logssb(y));
	if(y < S_logfsb.size())
	  total_logfsb(yall) = logspace_add2(total_logfsb(yall),S_logfsb(y));
	if(y < S_logCatch.size())
	  total_logCatch(yall) = logspace_add2(total_logCatch(yall),S_logCatch(y));
	if(y < S_logLand.size())
	  total_logLand(yall) = logspace_add2(total_logLand(yall),S_logLand(y));
	if(y < S_logtsb.size())
	  total_logtsb(yall) = logspace_add2(total_logtsb(yall),S_logtsb(y));
	if(y < S_logR.size())
	  total_logR(yall) = logspace_add2(total_logR(yall),S_logR(y));
	if(y < S_logfbar.size()){
	  total_nStocks(yall) = logspace_add2(total_nStocks(yall), Type(0.0));
	  total_logfbar(yall) = logspace_add2(total_logfbar(yall),S_logfbar(y));
	}
	if(y < S_logfbarL.size()){
	  // total_nStocks(yall) = logspace_add2(total_nStocks(yall), Type(0.0));
	  total_logfbarL(yall) = logspace_add2(total_logfbarL(yall),S_logfbarL(y));
	}
	for(int aall = 0; aall < maxAgeAll - minAgeAll + 1; ++aall){
	  int a = aall - (conf.minAge - minAgeAll);
	  if(a >= 0 && a < conf.maxAge - conf.minAge + 1 && y < S_logCatAge.cols()) 
	    total_logCatAge(aall,yall) = logspace_add2(total_logCatAge(aall,yall),S_logCatAge(a,y));
	}
      }
    }  
  }

  // fbar to averages
  total_logfbar -= total_nStocks;
  total_logfbarL -= total_nStocks;
  
  ADREPORT_F(total_logssb,of);
  ADREPORT_F(total_logfsb,of);
  ADREPORT_F(total_logCatch,of);
  ADREPORT_F(total_logLand,of);
  ADREPORT_F(total_logtsb,of);
  ADREPORT_F(total_logR,of);
  ADREPORT_F(total_logfbar,of);
  ADREPORT_F(total_logfbarL,of);
  ADREPORT_F(total_logCatAge,of);

  return;
}
