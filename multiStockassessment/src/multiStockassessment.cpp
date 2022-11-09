//  --------------------------------------------------------------------------
// Copyright (c) 2014, Anders Nielsen <an@aqua.dtu.dk>,    
// Casper Berg <cbe@aqua.dtu.dk>, Kasper Kristensen <kkr@aqua.dtu.dk>,
// Mollie Brooks <molbr@aqua.dtu.dk>,
// and Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>.
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//   * Redistributions of source code must retain the above copyright
//     notice, this list of conditions and the following disclaimer.
//   * Redistributions in binary form must reproduce the above copyright
//     notice, this list of conditions and the following disclaimer in the
//     documentation and/or other materials provided with the distribution.
//   * Neither the name of the assessment tool SAM nor the
//     names of its contributors may be used to endorse or promote products
//     derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL ANDERS NIELSEN, CASPER BERG OR KASPER 
// KRISTENSEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  --------------------------------------------------------------------------

#include "SAM.h"
#include "MSAM.h"
//"../inst/include/multiSAM.hpp"

template<class Type>
vector<Type> toLogProportion(vector<Type> x){
  vector<Type> r(x.size());
  r.setZero();
  Type lps = 0.0;
  for(int i = 0; i < x.size(); ++i){
    r(i) = x(i);
    lps = logspace_add2(lps,x(i));
  }
  r -= lps;
  return r;
}
  

// template<class VT, class Type>
// struct fake_data_indicator :
//   public data_indicator<VT, Type> {
//   fake_data_indicator() : data_indicator<VT, Type>(VT()){};
//   fake_data_indicator(VT obs) : data_indicator<VT, Type>(obs, true) {};
// };


extern "C" {


  SEXP vecpar2list(SEXP array){
    if(!Rf_isArray(array))Rf_error("Argument must be an array");
    return asSEXP(asCmoeVector<double>(array));
  }

  SEXP matpar2list(SEXP array){
    if(!Rf_isArray(array))Rf_error("Argument must be an array");
    return asSEXP(asCmoeMatrix<double>(array));
  }


}



template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRUCT(sam,sam_data);
  DATA_STRUCT(sharedObs,shared_obs);
  DATA_STRUCT(geneticsData, genetic_data);

  DATA_INTEGER(reportingLevel);
  
  ////////////////////////////////////////
  ////////// Multi SAM specific //////////
  ////////////////////////////////////////
  DATA_INTEGER(usePartialCor);
  DATA_INTEGER(maxYearAll);
  DATA_INTEGER(minYearAll);
  DATA_INTEGER(maxAgeAll);
  DATA_INTEGER(minAgeAll);
  DATA_STRUCT(cons, cov_constraints);
  DATA_MATRIX(X);
  DATA_MATRIX(XCon);

  DATA_INTEGER(ConConf);

  DATA_STRUCT(Xph, listMatrixFromR);

  // Related to residuals
  DATA_VECTOR(fake_obs);
  DATA_IVECTOR(fake_stock);
  DATA_IVECTOR(fake_indx);
  DATA_VECTOR_INDICATOR(fake_keep, fake_obs);
  DATA_INTEGER(doResiduals);

  // Create a keep indicator for each stock and for shared observations
  vector<data_indicator<vector<Type>,Type> > keep(sam.dataSets.size() + (sharedObs.hasSharedObs>0));
  
  for(int s = 0; s < sam.dataSets.size(); ++s){    
    data_indicator<vector<Type>,Type> keepTmp(sam.dataSets(s).logobs, true);
    keep(s) = keepTmp;
  }
  if(sharedObs.hasSharedObs){
    data_indicator<vector<Type>,Type> keepTmp(sharedObs.logobs, true);
    keep(keep.size()-1) = keepTmp;
  }
  
  if(doResiduals){
    for(int i = 0; i < fake_obs.size(); ++i){
      // Overwrite true observations with residual observations
      if(fake_stock(i) == sam.dataSets.size()){ // Shared observation
	sharedObs.logobs(fake_indx(i)) = fake_obs(i);
      }else{
	sam.dataSets(fake_stock(i)).logobs(fake_indx(i)) = fake_obs(i);
      }
      // Overwrite keep value
      keep(fake_stock(i))(fake_indx(i)) = fake_keep(i);
    }
  }

  PARAMETER_CMOE_VECTOR(logFpar);
  PARAMETER_CMOE_VECTOR(logQpow);
  PARAMETER_CMOE_VECTOR(logSdLogFsta);
  PARAMETER_CMOE_VECTOR(logSdLogN);
  PARAMETER_CMOE_VECTOR(logSdLogP);
  PARAMETER_CMOE_VECTOR(logSdLogObs);
  PARAMETER_CMOE_VECTOR(logSdLogTotalObs);
  PARAMETER_CMOE_VECTOR(transfIRARdist);
  PARAMETER_CMOE_VECTOR(sigmaObsParUS);
  PARAMETER_CMOE_VECTOR(rec_pars);
  PARAMETER_CMOE_VECTOR(itrans_rho);
  PARAMETER_CMOE_VECTOR(rhop);
  
  PARAMETER_CMOE_VECTOR(logScale);
  PARAMETER_CMOE_VECTOR(logitReleaseSurvival);
  PARAMETER_CMOE_VECTOR(logitRecapturePhi);
  PARAMETER_CMOE_VECTOR(logAlphaSCB);

  PARAMETER_CMOE_VECTOR(sepFalpha);
  PARAMETER_CMOE_VECTOR(sepFlogitRho);
  PARAMETER_CMOE_VECTOR(sepFlogSd);
  PARAMETER_CMOE_VECTOR(predVarObs);

  PARAMETER_CMOE_VECTOR(logPhiSW);
  PARAMETER_CMOE_VECTOR(logSdProcLogSW);
  PARAMETER_CMOE_VECTOR(meanLogSW);
  PARAMETER_CMOE_VECTOR(logSdLogSW);
  PARAMETER_CMOE_MATRIX(logPhiCW);
  PARAMETER_CMOE_VECTOR(logSdProcLogCW);
  PARAMETER_CMOE_VECTOR(meanLogCW);
  PARAMETER_CMOE_VECTOR(logSdLogCW);
  PARAMETER_CMOE_VECTOR(logPhiMO);
  PARAMETER_CMOE_VECTOR(logSdProcLogitMO);
  PARAMETER_CMOE_VECTOR(meanLogitMO);
  PARAMETER_CMOE_VECTOR(logSdMO);
  PARAMETER_CMOE_VECTOR(logPhiNM);
  PARAMETER_CMOE_VECTOR(logSdProcLogNM);
  PARAMETER_CMOE_VECTOR(meanLogNM);
  PARAMETER_CMOE_VECTOR(logSdLogNM);
  PARAMETER_CMOE_VECTOR(logXtraSd);
  PARAMETER_CMOE_VECTOR(initF);
  PARAMETER_CMOE_VECTOR(initN);

  // Forecast FMSY
  PARAMETER_VECTOR(logFScaleMSY);
  PARAMETER_VECTOR(implicitFunctionDelta);

  // YPR reference points
  // PARAMETER_VECTOR(logScaleFmsy);
  // PARAMETER_VECTOR(logScaleFmypyl);
  // PARAMETER_VECTOR(logScaleFmdy);
  // PARAMETER_VECTOR(logScaleFmax);
  // PARAMETER_VECTOR(logScaleF01);
  // PARAMETER_VECTOR(logScaleFcrash);
  // PARAMETER_VECTOR(logScaleFext);
  // PARAMETER_CMOE_VECTOR(logScaleFxPercent);
  // PARAMETER_VECTOR(logScaleFlim);
  // PARAMETER_CMOE_MATRIX(logScaleFmsyRange);

  PARAMETER_VECTOR(splinePenalty);


  PARAMETER_CMOE_MATRIX(logF); 
  PARAMETER_CMOE_MATRIX(logN);

  PARAMETER_CMOE_MATRIX(logSW);
  PARAMETER_CMOE_3DARRAY(logCW);  
  PARAMETER_CMOE_MATRIX(logitMO);
  PARAMETER_CMOE_MATRIX(logNM);

  PARAMETER_CMOE_MATRIX(logP);
 
  
  PARAMETER_CMOE_VECTOR(missing);
  
  ////////////////////////////////////////
  ////////// Multi SAM specific //////////
  ////////////////////////////////////////
  PARAMETER_VECTOR(RE);
  PARAMETER_VECTOR(betaCon); // (nages*nAreas) x (nages*nAreas)

  // PARAMETER_VECTOR(shared_logSdObs);
  PARAMETER_CMOE_VECTOR(shared_logFscale);
  // PARAMETER_VECTOR(shared_lfsMean);
  PARAMETER_VECTOR(shared_lfsSd);
  // PARAMETER_VECTOR(shared_lfsRho);
  PARAMETER_VECTOR(shared_logitMissingVulnerability);

  PARAMETER_CMOE_VECTOR(shared_phbeta);

  PARAMETER_CMOE_VECTOR(initLogN);
  PARAMETER_CMOE_VECTOR(initLogF);

  int nStocks = logF.cols();

  vector<paraSet<Type> > paraSets(nStocks);

  
  for(int s = 0; s < nStocks; ++s){
    paraSets(s).logFpar = logFpar.col(s);
    paraSets(s).logQpow = logQpow.col(s);  
    paraSets(s).logSdLogFsta = logSdLogFsta.col(s);
    paraSets(s).logSdLogN = logSdLogN.col(s);
    paraSets(s).logSdLogP = logSdLogP.col(s);
    paraSets(s).logSdLogObs = logSdLogObs.col(s);
    paraSets(s).logSdLogTotalObs = logSdLogTotalObs.col(s);
    paraSets(s).transfIRARdist = transfIRARdist.col(s);
    paraSets(s).sigmaObsParUS = sigmaObsParUS.col(s);
    paraSets(s).rec_pars = rec_pars.col(s);    
    paraSets(s).itrans_rho = itrans_rho.col(s);
    paraSets(s).rhop = rhop.col(s);
    
    paraSets(s).logScale = logScale.col(s);
    paraSets(s).logitReleaseSurvival = logitReleaseSurvival.col(s);
    paraSets(s).logitRecapturePhi = logitRecapturePhi.col(s);
    paraSets(s).logAlphaSCB = logAlphaSCB.col(s);
    
    paraSets(s).sepFalpha = sepFalpha.col(s);
    paraSets(s).sepFlogitRho = sepFlogitRho.col(s);
    paraSets(s).sepFlogSd = sepFlogSd.col(s);
    paraSets(s).predVarObs = predVarObs.col(s);
    // Forecast FMSY
    paraSets(s).logFScaleMSY = logFScaleMSY(s);
    paraSets(s).implicitFunctionDelta = implicitFunctionDelta(s);
    // Biopar    
    paraSets(s).logPhiSW = logPhiSW.col(s);
    paraSets(s).logSdProcLogSW = logSdProcLogSW.col(s);
    paraSets(s).meanLogSW = meanLogSW.col(s);
    paraSets(s).logSdLogSW = logSdLogSW.col(s);
    paraSets(s).logPhiCW = logPhiCW.col(s);
    paraSets(s).logSdProcLogCW = logSdProcLogCW.col(s);
    paraSets(s).meanLogCW = meanLogCW.col(s);
    paraSets(s).logSdLogCW = logSdLogCW.col(s);
    paraSets(s).logPhiMO = logPhiMO.col(s);
    paraSets(s).logSdProcLogitMO = logSdProcLogitMO.col(s);
    paraSets(s).meanLogitMO = meanLogitMO.col(s);
    paraSets(s).logSdMO = logSdMO.col(s);
    paraSets(s).logPhiNM = logPhiNM.col(s);
    paraSets(s).logSdProcLogNM = logSdProcLogNM.col(s);
    paraSets(s).meanLogNM = meanLogNM.col(s);
    paraSets(s).logSdLogNM = logSdLogNM.col(s);
    paraSets(s).logXtraSd = logXtraSd.col(s);
    paraSets(s).initF = initF.col(s);
    paraSets(s).initN = initN.col(s);
    if(splinePenalty.size() > 0)
      paraSets(s).splinePenalty = splinePenalty(s);
    // YPR reference points
    // paraSets(s).logScaleFmsy = logScaleFmsy(s);
    // paraSets(s).logScaleFmax = logScaleFmax(s);
    // paraSets(s).logScaleF01 = logScaleF01(s);
    // paraSets(s).logScaleFcrash = logScaleFcrash(s);
    // paraSets(s).logScaleFext = logScaleFext(s);
    // paraSets(s).logScaleFxPercent = logScaleFxPercent.col(s);
    // paraSets(s).logScaleFlim = logScaleFlim(s);
    // paraSets(s).logScaleFmsyRange = logScaleFmsyRange.col(s);
  }
  genetic_parameters<Type> genParSet;
  PARAMETER(gen_logKappaSpace); genParSet.logKappaSpace = gen_logKappaSpace;
  PARAMETER(gen_logKappaTime); genParSet.logKappaTime = gen_logKappaTime;
  PARAMETER_VECTOR(gen_corparST); genParSet.corparST = gen_corparST;
  PARAMETER_VECTOR(gen_logSdST); genParSet.logSdST = gen_logSdST;
  PARAMETER(gen_corparAge); genParSet.corparAge = gen_corparAge;
  PARAMETER_VECTOR(gen_corparTrip); genParSet.corparTrip = gen_corparTrip;
  PARAMETER_VECTOR(gen_logSdTrip); genParSet.logSdTrip = gen_logSdTrip;
  PARAMETER_ARRAY(gen_alleleFreq); genParSet.alleleFreq = gen_alleleFreq;
  PARAMETER_MATRIX(gen_dmScale); genParSet.dmScale = gen_dmScale;
 PARAMETER_MATRIX(gen_muLogP); genParSet.muLogP = gen_muLogP;

  PARAMETER_ARRAY(logGst);
  PARAMETER_MATRIX(logGtrip);


  // for(int s = 0; s < logF.cols(); ++s){  
  //   if(shared_logFscale.col(s).size() > 0){    
  //     matrix<Type> logF0 = logF.col(0);
  //     matrix<Type> logFs = logF.col(0);
  //     for(int i = 0; i < logF0.cols(); ++i){	  
  // 	for(int j = 0; j < logF0.rows(); ++j){	  
  // 	  //ans -= dnorm(logFs(j,i), logF0(j,i) + shared_logFscale.col(s)(i), Type(0.01), true);
  // 	  //ans -= dnorm(logFs(j,i), Type(0.0), Type(0.01), true);
  // 	  logFs(j,i) = logF0(j,i) + shared_logFscale.col(s)(i);
  // 	}
  //     }
  //     logF.col(s) = logFs;
  //   }
  // }


  // Spline penalty


  Type ans = 0.0; //negative log-likelihood

  ////////////////////////////////////
  /////////// Missing Obs ///////////
  //////////////////////////////////

 // patch missing 
  int idxmis = 0;
  for(int s = 0; s < nStocks; ++s)
    for(int i = 0; i < sam.dataSets(s).nobs; ++i){
      if(isNA(sam.dataSets(s).logobs(i))){
  	sam.dataSets(s).logobs(i) = missing(idxmis++);
      }    
    }

  // patch missing vulnerability keys
  idxmis = 0;
  for(int i = 0; i < sharedObs.keyFleetStock.rows(); ++i)
    for(int j = 0; j < sharedObs.keyFleetStock.cols(); ++j)
      if(isNA(sharedObs.keyFleetStock(i,j))){
  	sharedObs.keyFleetStock(i,j) = invlogit(shared_logitMissingVulnerability(idxmis++));
      }    
    
  
  // add wide prior for missing random effects, but _only_ when computing ooa residuals
  if(doResiduals){
    Type huge = 10.0;
    for (int i = 0; i < missing.size(); i++)
      ans -= dnorm(missing(i), Type(0.0), huge, true);
  } 


  
  ////////////////////////////////////
  /////////// Recruitment ///////////
  //////////////////////////////////
  vector<Recruitment<Type> > recruits(nStocks);
  for(int s = 0; s < nStocks; ++s){
    recruits(s) = makeRecruitmentFunction(sam.confSets(s), paraSets(s));
  }

  ////////////////////////////////////////
  /////////// Prepare forecast //////////
  //////////////////////////////////////

  for(int s = 0; s < nStocks; ++s){
    // Calculate forecast
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    // Resize arrays
    prepareForForecast(sam.forecastSets(s), sam.dataSets(s), sam.confSets(s), paraSets(s), logFa, logNa, recruits(s));
  }
     
  ofall<Type> ofAll(nStocks);

  /////////////////////////////////////
  ////////// Spline Penalty //////////
  ///////////////////////////////////

  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    // array<Type> logNa = getArray(logN, s);
    // array<Type> logFa = getArray(logF, s);
    // data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);
    ans += nllSplinePenalty(ds, cs, ps, &of);
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }

  //////////////////////////////////
  ////////// Bio PROCESS //////////
  ////////////////////////////////

  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    array<Type> logSWa = getArray(logSW, s);
    array<Type> logCWa = getArray(logCW, s);
    array<Type> logitMOa = getArray(logitMO, s);
    array<Type> logNMa = getArray(logNM, s);
   
    // data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);

    ans += nllSW(logSWa, ds, cs, ps, &of);
    ans += nllCW(logCWa, ds, cs, ps, &of);
    ans += nllMO(logitMOa, ds, cs, ps, &of);
    ans += nllNM(logNMa, ds, cs, ps, &of);      

    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }

  

  //////////////////////////////////
  /////////// Components //////////
  ////////////////////////////////

  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    array<Type> logPa = getArray(logP, s);
    // Resize arrays
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);
    data_indicator<vector<Type>,Type> keepTmp = keep(s);
    ans += nllP(cs, ps, logPa, keepTmp, &of);
  }
     
  
  /////////////////////////////////////////
  ////////// F PRE-CALCULATIONS //////////
  ///////////////////////////////////////

  vector<bool> hasPH(Xph.size());
  hasPH.setConstant(false);
  vector<matrix<Type> > phPred(Xph.size());
  for(int s = 0; s < phPred.size(); ++s){
    if(shared_phbeta.col(s).size() > 0){
      vector<Type> v0 = Xph(s) * shared_phbeta.col(s);
      if(v0.size() != (maxAgeAll - minAgeAll + 1)*(maxYearAll - minYearAll + 1))
	Rf_error("Wrong size in proportional hazard");
      phPred(s) = asMatrix(v0, maxAgeAll - minAgeAll + 1, maxYearAll - minYearAll + 1); // Age x stock
      hasPH(s) = true;
    }else{
      matrix<Type> m0(0,0);
      m0.setZero();
      phPred(s) = m0;
    }
  }
  REPORT(phPred);
    
  // Initial parameter contribution
  for(int s = 0; s < nStocks; ++s){
    if(initLogF.col(s).size() > 0){
      if(shared_logFscale.col(s).size() == 0){
	for(int a = 0; a < initLogF.col(s).size(); ++a)
	  ans -= dnorm(logF.col(s)(a,0), initLogF.col(s)(a), Type(0.01), true);
      }else{
	ans -= dnorm(shared_logFscale.col(s)(0), initLogF.col(s)(0), Type(0.01), true);
      }
    }
  }

  // bool overwriteF = false;
  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);
    if(shared_logFscale.col(s).size() == 0 && !hasPH(s)){ // Not using shared logF selectivity or proportional hazard
      // SAVE FOR LATER
      // ans += nllF(ds, cs, ps, sam.forecastSets(s), logFa, keepTmp, &of);
      // ofAll.addToReport(of.report,s);
      // moveADREPORT(&of,this,s);
      // // If simulate -> move grab new logF values and move them to the right place!         
      // logF.col(s) = logFa.matrix();
    }else{			// Using shared logF selectivity
      matrix<Type> logF0 = logF.col(0);
      matrix<Type> logFs = logF.col(s);      
      for(int i = 0; i < logFs.cols(); ++i){ // Loop over time
	Type slfs = 0.0;
	if(shared_logFscale.col(s).size() > 0){
	  if(i > 0){
	    //ans -= dnorm(shared_logFscale.col(s)(i) - shared_lfsMean(s-1), invlogit(shared_lfsRho(s-1)) * (shared_logFscale.col(s)(i-1) - shared_lfsMean(s-1)), exp(shared_lfsSd(s-1)), true);
	    ans -= dnorm(shared_logFscale.col(s)(i), shared_logFscale.col(s)(i-1), exp(shared_lfsSd(s-1)), true);
	  }else if(initLogF.col(s).size() == 0){
	    //Type logsd = shared_lfsSd(s-1) - 0.5 * logspace_add2(log(2.0)+shared_lfsRho(s-1),Type(0.0)) - logspace_add2(shared_lfsRho(s-1),Type(0.0));
	    //ans -= dnorm(shared_logFscale.col(s)(0),shared_lfsMean(s-1),exp(logsd),true);
	    ans -= dnorm(shared_logFscale.col(s)(0),Type(0.0),Type(2.0),true);
	  }
	  slfs = shared_logFscale.col(s)(i);
	}
	for(int j = 0; j < logF0.rows(); ++j){	  
	  Type lfPred = logF0(j,i);
	  if(hasPH(s)){
	    Type phV = 0.0;
	    Type phD = 0.0;
	    for(int a = 0; a < phPred.rows(); ++a){
	      if(cs.keyLogFsta(0,a) == j){
		phV += phPred(s)(a,i);
		phD += 1.0;
	      }
	    }
	    if(phD > 0)
	      lfPred += phV / phD;
	  }
	  //ans -= dnorm(logFs(j,i), lfPred + slfs, Type(0.01), true);
	  ans -= dnorm(logFs(j,i), Type(0.0), Type(1.0), true);
	  logFs(j,i) = lfPred + slfs;
	}
	// Implement simulation of logF with shared selectivity!
      }
      logF.col(s) = logFs;
      //overwriteF = true;
    }
  }

  // if(overwriteF){
  for(int s = 0; s < logF.cols(); ++s){
    oftmp<Type> of(this->do_simulate);
    matrix<Type> logFs = logF.col(s);
    ADREPORT_F(logFs,(&of));
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }
  // }


 
  ////////////////////////////////////////
  /////////// Mortality /////////////////
  //////////////////////////////////////
  vector<MortalitySet<Type> > mortalities(nStocks);
  for(int s = 0; s < nStocks; ++s){
    array<Type> logFa = getArray(logF, s);
    mortalities(s) = MortalitySet<Type>(sam.dataSets(s), sam.confSets(s), paraSets(s), logFa);
  }

  ////////////////////////////////////////
  /////////// Calculate forecast ////////
  //////////////////////////////////////

  for(int s = 0; s < nStocks; ++s){
    // Calculate forecast
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    // Resize arrays
    sam.forecastSets(s).calculateForecast(logFa,logNa, sam.dataSets(s), sam.confSets(s), paraSets(s), recruits(s), mortalities(s));
  }

  ////////////////////////////////
  ////////// F PROCESS //////////
  //////////////////////////////


   for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);
    if(shared_logFscale.col(s).size() == 0 && !hasPH(s)){ // Not using shared logF selectivity or proportional hazard
      ans += nllF(ds, cs, ps, sam.forecastSets(s), logFa, keepTmp, &of);
      ofAll.addToReport(of.report,s);
      moveADREPORT(&of,this,s);
      // If simulate -> move grab new logF values and move them to the right place!         
      logF.col(s) = logFa.matrix();
    }
   }

  
  ////////////////////////////////
  ////////// N PROCESS //////////
  //////////////////////////////

  // Initial parameter contribution
  for(int s = 0; s < nStocks; ++s){
    if(initLogN.col(s).size() == logN.col(s).rows()){ // Initial value per age
      for(int a = 0; a < initLogN.col(s).size(); ++a)
	ans -= dnorm(logN.col(s)(a,0), initLogN.col(s)(a), Type(0.01), true);   
    }else if(initLogN.col(s).size() == 1){			// Initial value for "back-in-time" recruitment assuming first year M and F were constant back in time
      //predN(0) = initLogN.col(s)(0);
      confSet conf = sam.confSets(s);
      dataSet<Type> dat = sam.dataSets(s);
      paraSet<Type> par = paraSets(s);
      array<Type> logNa = getArray(logN, s);
      array<Type> logFa = getArray(logF, s);

      // Type logssb0 = ssbi(dat, conf, logNa, logFa, 0, true);
      // ans -= dnorm(logssb0, initLogN.col(s)(0), Type(1.0), true);
      ans -= dnorm(logNa(0,0), initLogN.col(s)(0), Type(0.01), true);

      for(int a = 1; a < logNa.rows(); ++a){
      	Type m = logNa(a-1,0) - dat.natMor(0,a-1);
      	if(conf.keyLogFsta(0,a-1)>(-1))
      	  m -= exp(logFa(conf.keyLogFsta(0,a-1),0));	
      	// ans -= dnorm((logNa(a,0)), (m), exp(par.logSdLogN(conf.keyVarLogN(a))), true); //
      	ans -= dnorm((logNa(a,0)), (m), Type(0.01), true); //
      }
    }
  }      
  //////////////////////////////////////////

  int nAreas = sam.dataSets.size();

  int nages = (maxAgeAll - minAgeAll + 1);

  matrix<Type> AlphaCon(nages*nAreas,nages*nAreas);
  AlphaCon.setZero();
  vector<Type> betaConX = (matrix<Type>)XCon * (vector<Type>)betaCon;
  int k = 0;
  if(betaConX.size() > 0){
    for(int i = 0; i < AlphaCon.cols(); ++i){
      for(int j = 0; j < AlphaCon.rows(); ++j){
	if(ConConf == 1 && i == j){
	  AlphaCon(j,i) = 0.0;
	}else if(ConConf == 2 && i == j){
	  AlphaCon(j,i) = 1.0;
	}else{
	  AlphaCon(j,i) = betaConX(k++);
	}
      }
    }
  }
  REPORT(AlphaCon);
  REPORT(betaConX);
  
  matrix<Type> A(nages*nAreas,nages*nAreas);
  A.setZero();
  // Create matrix A of variances to scale

  for(int i = 0; i < A.rows(); ++i){ // Loop over indices
    // Area
    int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
    // Age index = age - minAge
    int a = i % nages;
    if(a + minAgeAll >= sam.confSets(s).minAge && a + minAgeAll <= sam.confSets(s).maxAge){
      Type lsdln = paraSets(s).logSdLogN(sam.confSets(s).keyVarLogN(a - (sam.confSets(s).minAge - minAgeAll)));
      A(i,i) = exp(lsdln);
    }else{
      A(i,i) = Type(1.0);
    }
  }
  REPORT(A);

  matrix<Type> L = constructL(nages,nAreas,RE, X,cons);
    
  // Get covariance matrix with arbitrary scale
  matrix<Type> Sigma = L * L.transpose();
  //Sigma = cov2cor(Sigma);
  matrix<Type> SigmaTmp = Sigma;
    
  if(usePartialCor){
    Sigma = atomic::matinv(Sigma); //cov2cor((matrix<Type>)Sigma.inverse()); //pcor2cor(Sigma);
  }

  Sigma = cov2cor(Sigma);
  
  matrix<Type> ncov = (matrix<Type>)(A * Sigma) * A;
  REPORT(Sigma);
  REPORT(SigmaTmp);
  REPORT(L);
  REPORT(ncov)

    ////////////////////////////////////
    ////// CALCULATE CONTRIBUTION //////
    ////////////////////////////////////

    // add wide prior for first state, but _only_ when computing ooa residuals
    if(doResiduals){
      Type huge = 10.0;
      for(int s = 0; s < nAreas; ++s){
	for (int i = 0; i < logN.col(s).rows(); i++){
	  ans -= dnorm((logN.col(s))(i, 0), Type(0.0), huge, true);
	}
      }
    } 
    
  density::MVNORM_t<Type> neg_log_densityN(ncov);
  // Loop over time
  for(int yall = 0; yall < maxYearAll - minYearAll + 1; ++yall){
    // Handle simulation of F for forecast and HCR
    SIMULATE{
      // Simulate new F for forecast and HCR
     for(int s = 0; s < nAreas; ++s){
       array<Type> logNa = getArray(logN, s);
       array<Type> logFa = getArray(logF, s);
       forecastSet<Type> fc = sam.forecastSets(s);
       if(fc.nYears > 0){
	 dataSet<Type> ds = sam.dataSets(s);
	 confSet cs = sam.confSets(s);
	 paraSet<Type> ps = paraSets(s);
	 int y = yall - CppAD::Integer(ds.years(0) - minYearAll);
	 if(y > 0 && y < ds.noYears + fc.nYears){
	   matrix<Type> fvar = get_fvar(ds, cs, ps, logFa);
	   MVMIX_t<Type> neg_log_densityF(fvar,Type(cs.fracMixF));
	   int nYears = fc.nYears;
	   int fi = y - fc.forecastYear.size() + nYears;
	   // Update forecast
	   if(fi > 0){
	     fc.updateForecast(fi, logFa, logNa, ds, cs, ps, recruits(s), mortalities(s));
	     // Simulate F
	     // int forecastIndex = CppAD::Integer(dat.forecast.forecastYear(i))-1;
	     if(fc.simFlag(0) == 0){
	       Type timeScale = fc.forecastCalculatedLogSdCorrection(fi);
	       logF.col(s).col(y) = (vector<Type>)fc.forecastCalculatedMedian.col(fi) + neg_log_densityF.simulate() * timeScale;	 
	     }	   
	   }	 
	 }
       }
     }
    }
    
    // Vector for predictions
    vector<Type> lastN(ncov.rows());
    lastN.setZero();
    vector<Type> predN(ncov.rows());
    predN.setZero();
    vector<Type> newN(ncov.rows());
    newN.setZero();
    // Determine which ages to use with keep vector
    vector<Type> keepN(ncov.rows());
    keepN.setZero();
    // Loop over stocks
    for(int s = 0; s < nAreas; ++s){
      dataSet<Type> ds = sam.dataSets(s);
      confSet cs = sam.confSets(s);
      paraSet<Type> ps = paraSets(s);
      int ageOffset = sam.confSets(s).minAge - minAgeAll;
      int y = yall - CppAD::Integer(ds.years(0) - minYearAll);
      
      array<Type> logNa = getArray(logN, s);
      array<Type> logFa = getArray(logF, s);
      if(y > 0 && y < ds.noYears + sam.forecastSets(s).nYears){
	vector<Type> predNnz = predNFun(ds, cs, ps, logNa, logFa, recruits(s), mortalities(s), (int)y);
	keepN.segment(s * nages + ageOffset,predNnz.size()) = 1.0;
	predN.segment(s * nages + ageOffset,predNnz.size()) = predNnz;
	newN.segment(s * nages + ageOffset,predNnz.size()) = logNa.col(y);
	lastN.segment(s * nages + ageOffset,predNnz.size()) = logNa.col(y-1);	
      }
    }
    predN -= AlphaCon * lastN;
    if(keepN.sum()>0){
      // forecast correction to recruitment
      vector<Type> Nscale(predN.size());
      Nscale.setZero();
      Nscale += 1.0;
      vector<Type> predNTmp = predN;

      for(int s = 0; s < nAreas; ++s){
	dataSet<Type> ds = sam.dataSets(s);
	int ageOffset = sam.confSets(s).minAge - minAgeAll;
	int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	if(y > 0 &&
	   y < ds.noYears + sam.forecastSets(s).nYears &&
	   sam.forecastSets(s).nYears > 0 &&
	   sam.forecastSets(s).recModel(CppAD::Integer(sam.forecastSets(s).forecastYear(y))-1) != sam.forecastSets(s).asRecModel &&
	   sam.forecastSets(s).forecastYear(y) > 0){
	  Nscale(s * nages + ageOffset) = sqrt(sam.forecastSets(s).logRecruitmentVar) / sqrt(ncov(s * nages + ageOffset,s * nages + ageOffset));
	  predNTmp(s * nages + ageOffset) = sam.forecastSets(s).logRecruitmentMedian;
	}
      }
      ans+=neg_log_densityN((newN-predNTmp) / Nscale, keepN) + (keepN * log(Nscale)).sum();	  
    }// else{	// end forecast correction to recruitment
    //   ans+=neg_log_densityN(newN-predN, keepN);
    // }
      
    SIMULATE{
      /* Plan:
	 1) If all sam.confSets(s).simFlag(1) == 1, skip the simulation
	 2) Get conditional distribution of stocks with simFlag == 0 (and keepN == 1)
	 3) Simulate from marginal distribution
	 4) Insert into logN at the right places
      */
      // 1) Check if any simFlags are 0
      bool doSim = false;
      for(int s = 0; s < nAreas; ++s){
	int nYears = sam.forecastSets(s).nYears;
	if(nYears > 0){
	  int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	  int fi = y - sam.forecastSets(s).forecastYear.size() + nYears;
	  if(sam.confSets(s).simFlag(1) == 0 ||
	     (fi > 0 && sam.forecastSets(s).simFlag(1) == 0)){
	    doSim = true;
	    break;
	  }
	}else{
	  if(sam.confSets(s).simFlag(1) == 0){
	    doSim = true;
	    break;
	  }
	}
      }
      if(doSim){
	// 2) Get conditional distribution of stocks with simFlag == 0 (and keepN = 1)
	vector<int> notCondOn(ncov.cols());
	notCondOn.setZero();
	for(int i = 0; i < notCondOn.size(); ++i){
	  int s = (int)i / (int)nages;
	  int nYears = sam.forecastSets(s).nYears;
	  if(nYears > 0){
	    int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	    int fi = y - sam.forecastSets(s).forecastYear.size() + nYears;
	    if((sam.confSets(s).simFlag(1) == 0 ||
		(fi > 0 && sam.forecastSets(s).simFlag(1) == 0)) &&
	       keepN(i) == 1){
	      notCondOn(i) = 1;
	    }
	  }else if(sam.confSets(s).simFlag(1) == 0 && keepN(i) == 1){
	    notCondOn(i) = 1;
	  }
	}

	int nAll = ncov.cols();
	int nNotCond = notCondOn.sum();
	int nCond = nAll - nNotCond;
	if(nNotCond > 0){ // Only do this if there are any not conditional
	  vector<int> ccond(nNotCond);
	  vector<int> cond(nCond);
	  int k1 = 0; int k2 = 0;
	  for(int i = 0; i < notCondOn.size(); ++i){
	    if(notCondOn(i) == 0){
	      cond(k1++) = i;
	    }else{
	      ccond(k2++) = i;
	    }
	  }
	  matrix<Type> Sigma_NN(nNotCond,nNotCond);
	  matrix<Type> Sigma_CC(nCond,nCond);
	  matrix<Type> Sigma_NC(nNotCond,nCond);
	  matrix<Type> Sigma_CN(nCond,nNotCond);

	  for(int i = 0; i < ccond.size(); ++i)
	    for(int j = 0; j < ccond.size(); ++j)
	      Sigma_NN(i,j) = ncov(ccond(i),ccond(j));

	  for(int i = 0; i < cond.size(); ++i)
	    for(int j = 0; j < cond.size(); ++j)
	      Sigma_CC(i,j) = ncov(cond(i),cond(j));
	  matrix<Type> Sigma_CC_inv = Sigma_CC.inverse();

	  for(int i = 0; i < ccond.size(); ++i)
	    for(int j = 0; j < cond.size(); ++j)
	      Sigma_NC(i,j) = ncov(ccond(i),cond(j));
	  
	  matrix<Type> meanCorrectionMat = Sigma_NC * Sigma_CC_inv;

	  for(int i = 0; i < cond.size(); ++i)
	    for(int j = 0; j < ccond.size(); ++j)
	      Sigma_CN(i,j) = ncov(cond(i),ccond(j));

	  matrix<Type> newSigma = Sigma_NN - (matrix<Type>)(meanCorrectionMat * Sigma_CN);

	  // 3) Simulate marginal distribution
	  vector<Type> avec(nCond);
	  vector<Type> muNew(nNotCond);
	  k1 = 0; k2 = 0;
	  for(int i = 0; i < notCondOn.size(); ++i){
	    if(notCondOn(i) == 0){
	      if(keepN(i) == 0){
		avec(k1++) = 0.0;
	      }else{
		int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
		// Age index = age - minAge
		int a = i % nages;
		int ageOffset = sam.confSets(s).minAge - minAgeAll;
		int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);	
		avec(k1++) = logN.col(s)(a - ageOffset,y) - predN(i);
	      }
	    }else{
	      muNew(k2++) = predN(i);
	    }
	  }
	  vector<Type> simRes = muNew + (vector<Type>)(meanCorrectionMat * avec) + density::MVNORM(newSigma).simulate();
	  
	  // 4) Insert into logN at the right places
	  k1 = 0;
	  for(int i = 0; i < notCondOn.size(); ++i){
	    if(notCondOn(i) == 1){
	      int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
	      // Age index = age - minAge
	      int a = i % nages;
	      int ageOffset = sam.confSets(s).minAge - minAgeAll;
	      int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	      logN.col(s)(a - ageOffset,y) = simRes(k1++);
	    }
	  }
	}
      }// End doSim
    } // End simulate
  } // End year loop

  ///////////////////////////////////
  ////////// OBSERVATIONS //////////
  /////////////////////////////////

  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    array<Type> logPa = getArray(logP, s);
    data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);
    ans += nllObs(ds, cs, ps, sam.forecastSets(s), logNa, logFa, logPa, recruits(s), mortalities(s), keepTmp, reportingLevel,  &of);
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }

  // ans += sharedObservation(sharedObs,
  // 			   sam.dataSets,
  // 			   sam.confSets,
  // 			   paraSets,
  // 			   sam.forecastSets,
  // 			   logF,
  // 			   logN,
  // 			   logP,
  // 			   mortalities,
  // 			   // shared_logSdObs,
  // 			   minYearAll,
  // 			   minAgeAll,
  // 			   keep(keep.size()-1),
  // 			   this);

  ///////////////////////////////
  ////////// GENETICS //////////
  /////////////////////////////

  // ans += nllGenetics(sharedObs,
  // 		     sam.dataSets,
  // 		     sam.confSets,
  // 		     paraSets,
  // 		     sam.forecastSets,
  // 		     logF,
  // 		     logN,
  // 		     logP,
  // 		     mortalities,
  // 		     genParSet,
  // 		     geneticsData,
  // 		     logGst,
  // 		     logGtrip,
  // 		     maxAgeAll,
  // 		     minAgeAll);

  ///////////////////////////////////////
  ////////// REFERENCE POINTS //////////
  /////////////////////////////////////

  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
    logNa = logN.col(s);
    array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
    logFa = logF.col(s);
    data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);
    reportDeterministicReferencePoints(ds, cs, ps, logNa, logFa, recruits(s), sam.referencepointLists(s), &of);
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }

  ////////////////////////////////////////////////////////
  ////// ADD REPORTED OBJECTS FROM stockassessment //////
  //////////////////////////////////////////////////////

  if(isDouble<Type>::value && this->current_parallel_region<0)
    moveREPORT(this->report,ofAll.report);



  getTotals(sam.dataSets,
	    sam.confSets,
	    paraSets,
	    logF,
	    logN,
	    mortalities,
	    minYearAll,
	    maxYearAll,
	    minAgeAll,
	    maxAgeAll,
	    this);

  return ans;
}
