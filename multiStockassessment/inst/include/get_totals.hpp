
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
  vector<Type> total_logfbarL(nYear);
  total_logfbarL.setConstant(R_NegInf);

  vector<Type> total_nStocks(nYear);
  total_nStocks.setConstant(R_NegInf);
  
    
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
	if(y < S_logfbarL.size())
	  total_logfbarL(yall) = logspace_add2(total_logfbarL(yall),S_logfbarL(y));
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
