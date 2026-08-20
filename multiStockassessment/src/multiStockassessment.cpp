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


// template<class Type>
// tmbutils::array<Type> toArray(matrix<Type> x){
//   vector<int> ardi(2);
//   ardi << x.rows(), x.cols();
//   vector<Type> a = x.vec();
//   array<Type> y(a,ardi);
//   return y;
// }

template<class Type>
tmbutils::array<Type> toArray(matrix<Type> x){
  array<Type> y(x.rows(), x.cols());
  y = x;
  return y;
}





// template<class Type>
// void reportMort(MortalitySet<Type> mort, objective_function<Type>* of){
//   REPORT_F(mort.cumulativeHazard,of);
//   REPORT_F(mort.cumulativeHazard_F,of);
//   REPORT_F(mort.logFleetSurvival_before,of);
//   REPORT_F(mort.fleetCumulativeIncidence,of);
//   REPORT_F(mort.otherCumulativeIncidence,of);
//   REPORT_F(mort.ssbSurvival_before,of);
//   REPORT_F(mort.Fseason,of);
//   vector<Type> brk(mort.activeHazard_breakpoints);
//   REPORT_F(brk,of);
//   REPORT_F(mort.activeHazard_season,of);
//   REPORT_F(mort.activeHazard_F,of);
//   REPORT_F(mort.activeHazard_M,of);
//   REPORT_F(mort.Hazard_breakpoints,of);
//   REPORT_F(mort.CIF_F_breakpoints,of);
//   REPORT_F(mort.CIF_M_breakpoints,of);
//   return;
// }

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
  // SIMULATE{
  //   GetRNGstate();
  // }
  DATA_STRUCT(sam,sam_data);
  DATA_STRUCT(sharedObs,shared_obs);
  DATA_STRUCT(geneticsData, genetic_data);

  DATA_INTEGER(reportingLevel);
  // 0: Not for mohn's rho (i.e., multi stock assessment)  
  // 1: independence for simplicity
  // 2: quadratic form correction
  // 3: replace data with random effects
  DATA_INTEGER(mohn);		

  
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

  DATA_STRUCT(Xph, listArrayFromR);
  /*
    -1: Same F for all stocks
    0: Not shared, a RW by stock
    1: RW for first stock, others used the first scaled by vector AR1
    2: RW for first stock, others used the first scaled by scalar RW
    3: RW for first stock, others are scaled by log-linear relation
    4: RW for first stock, others are scaled by log-linear relation + scalar RW
   */
  DATA_INTEGER(shared_F_type);
  DATA_INTEGER(shared_Fseason_type);
  DATA_INTEGER(skip_stock_observations);

  DATA_IMATRIX(stockAreas); // areas x stocks. 1 if active, 0 otherwise
  DATA_IMATRIX(keyContamination);
  
  // Related to residuals
  DATA_VECTOR(fake_obs);
  DATA_IVECTOR(fake_stock);
  DATA_IVECTOR(fake_indx);
  DATA_VECTOR_INDICATOR(fake_keep, fake_obs);
  DATA_INTEGER(fake_includeOthers);
  DATA_INTEGER(doResiduals);
  DATA_INTEGER(useFakeObs);

  // Create a keep indicator for each stock and for shared observations
  vector<data_indicator<vector<Type>,Type> > keep(sam.dataSets.size() + (sharedObs.hasSharedObs>0));
  
  for(int s = 0; s < sam.dataSets.size(); ++s){    
    data_indicator<vector<Type>,Type> keepTmp(sam.dataSets(s).logobs, true);
    if(fake_includeOthers == 0)
      keepTmp.setZero();
    keepTmp.osa_flag = (bool)doResiduals;
    keep(s) = keepTmp;
  }
  if(sharedObs.hasSharedObs){
    data_indicator<vector<Type>,Type> keepTmp(sharedObs.logobs, true);
    if(fake_includeOthers == 0)
      keepTmp.setZero();
    keepTmp.osa_flag = (bool)doResiduals;
    keep(keep.size()-1) = keepTmp;
  }
  
  if(doResiduals || useFakeObs){
    for(int i = 0; i < fake_obs.size(); ++i){
      // Overwrite true observations with residual observations
      if(sharedObs.hasSharedObs && (fake_stock(i) == sam.dataSets.size())){ // Shared observation
	sharedObs.logobs(fake_indx(i)) = fake_obs(i);
      }else{
	sam.dataSets(fake_stock(i)).logobs(fake_indx(i)) = fake_obs(i);
      }
      // Overwrite keep value
      if(doResiduals)
	keep(fake_stock(i))(fake_indx(i)) = fake_keep(i);
    }
  }

    
  
  PARAMETER_CMOE_VECTOR(logFpar);
  PARAMETER_CMOE_VECTOR(logQpow);
  PARAMETER_CMOE_VECTOR(logSdLogFsta);
  PARAMETER_CMOE_VECTOR(muF);
  PARAMETER_CMOE_VECTOR(trans_rho_F);
  PARAMETER_CMOE_VECTOR(boundF_kappa);
  PARAMETER_CMOE_VECTOR(boundF_alpha);
  PARAMETER_CMOE_VECTOR(boundF_tau);
  PARAMETER_CMOE_VECTOR(boundFTAC_kappa)
  PARAMETER_CMOE_VECTOR(boundFTAC_alpha)
  PARAMETER_CMOE_VECTOR(logSdLogN);
  PARAMETER_CMOE_VECTOR(logSdLogP);
  PARAMETER_CMOE_VECTOR(logSdLogObs);
  PARAMETER_CMOE_VECTOR(logSdLogTotalObs);
  PARAMETER_CMOE_VECTOR(transfIRARdist);
  PARAMETER_CMOE_VECTOR(sigmaObsParUS);
  PARAMETER_CMOE_VECTOR(rec_pars);
  PARAMETER_CMOE_VECTOR(rec_transphi);
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
  PARAMETER_CMOE_VECTOR(recVarScalePar);
  PARAMETER_VECTOR(logFecundityScaling);
  PARAMETER_CMOE_VECTOR(logSpawningQuality);
  PARAMETER_CMOE_VECTOR(Wsigma);

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
  PARAMETER_CMOE_MATRIX(Mbeta);
  PARAMETER_CMOE_VECTOR(logSdLogNM);
  PARAMETER_CMOE_VECTOR(logXtraSd);
  PARAMETER_CMOE_VECTOR(initF);
  PARAMETER_CMOE_VECTOR(initN);
  PARAMETER_CMOE_VECTOR(scaleMpars);
  PARAMETER_CMOE_VECTOR(cp_m);
  PARAMETER_CMOE_VECTOR(cp_logk);
  PARAMETER_CMOE_VECTOR(cp_loga);
  PARAMETER_CMOE_VECTOR(cp_logb);


  PARAMETER_CMOE_MATRIX(seasonMu);
  PARAMETER_CMOE_VECTOR(seasonLogitRho);
  PARAMETER_CMOE_VECTOR(seasonLogSd);

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

  PARAMETER_CMOE_3DARRAY(logitFseason);
  
  PARAMETER_CMOE_VECTOR(missing);


   
  ////////////////////////////////////////
  ////////// Multi SAM specific //////////
  ////////////////////////////////////////

  
  PARAMETER_VECTOR(RE);
  PARAMETER_VECTOR(betaCon); // (nages*nAreas) x (nages*nAreas)

  // PARAMETER_VECTOR(shared_logSdObs);  
  PARAMETER_CMOE_MATRIX(shared_logFscale);
  PARAMETER_MATRIX(shared_lfsMean);
  PARAMETER_VECTOR(shared_lfsSd);
  PARAMETER_VECTOR(shared_lfsRho);
  PARAMETER_VECTOR(shared_logitMissingVulnerability);
  PARAMETER_VECTOR(shared_missingObs);

  PARAMETER_CMOE_VECTOR(shared_phbeta);
  
  PARAMETER_CMOE_VECTOR(initLogN);
  PARAMETER_CMOE_VECTOR(initLogF);

  PARAMETER_VECTOR(mohn_fake_obs);

  int nStocks = logF.cols();

  vector<paraSet<Type> > paraSets(nStocks);


  
  for(int s = 0; s < nStocks; ++s){
    paraSets(s).logFpar = logFpar.col(s);
    paraSets(s).logQpow = logQpow.col(s);  
    paraSets(s).logSdLogFsta = logSdLogFsta.col(s);
    paraSets(s).muF = muF.col(s);
    paraSets(s).trans_rho_F = trans_rho_F.col(s);
    paraSets(s).boundF_kappa = boundF_kappa.col(s);
    paraSets(s).boundF_alpha = boundF_alpha.col(s);
    paraSets(s).boundFTAC_kappa = boundFTAC_kappa.col(s);
    paraSets(s).boundFTAC_alpha = boundFTAC_alpha.col(s);
    paraSets(s).boundF_tau = boundF_tau.col(s);
    paraSets(s).logSdLogN = logSdLogN.col(s);
    paraSets(s).logSdLogP = logSdLogP.col(s);
    paraSets(s).logSdLogObs = logSdLogObs.col(s);
    paraSets(s).logSdLogTotalObs = logSdLogTotalObs.col(s);
    paraSets(s).transfIRARdist = transfIRARdist.col(s);
    paraSets(s).sigmaObsParUS = sigmaObsParUS.col(s);
    paraSets(s).rec_pars = rec_pars.col(s);
    paraSets(s).rec_transphi = rec_transphi.col(s);
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
    paraSets(s).recVarScalePar = recVarScalePar.col(s);
    paraSets(s).logFecundityScaling = logFecundityScaling(s);
    paraSets(s).logSpawningQuality = logSpawningQuality.col(s);
    paraSets(s).Wsigma = Wsigma.col(s);
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
    paraSets(s).Mbeta = Mbeta.col(s);
    paraSets(s).logSdLogNM = logSdLogNM.col(s);
    paraSets(s).logXtraSd = logXtraSd.col(s);
    paraSets(s).initF = initF.col(s);
    paraSets(s).initN = initN.col(s);
    paraSets(s).scaleMpars = scaleMpars.col(s);
    paraSets(s).cp_m = cp_m.col(s);
    paraSets(s).cp_logk = cp_logk.col(s);
    paraSets(s).cp_loga = cp_loga.col(s);
    paraSets(s).cp_logb = cp_logb.col(s);
    if(splinePenalty.size() > 0)
      paraSets(s).splinePenalty = splinePenalty(s);
    paraSets(s).seasonMu = seasonMu.col(s);
    paraSets(s).seasonLogitRho = seasonLogitRho.col(s);
    paraSets(s).seasonLogSd = seasonLogSd.col(s);
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
  PARAMETER_CMOE_MATRIX(gen_logitConfusionMatrix);
  // Convert from (identifiable) logit to (full) log
  if(gen_logitConfusionMatrix.size() > 0){
    // Rcout << "There is a confusion matrix!\n";
    // Rcout << gen_logitConfusionMatrix.cols() << "\n";
    vector<matrix<Type>> logCM(gen_logitConfusionMatrix.cols());
    // Rcout << logCM.size() << "\n";
    for(int i = 0; i < logCM.size(); ++i){
      // Rcout << "\t" << i << "\n";
      // a column per genetic stock, a row per possible classification. Rows must sum to 1
      matrix<Type> v0 = gen_logitConfusionMatrix.col(i);
      // Rcout << "\t" << v0.rows() << " x " << v0.cols() << "\n";
      matrix<Type> v1(v0.rows() + 1, v0.cols());
      // Rcout << "\t" << v1.rows() << " x " << v1.cols() << "\n";
      v1.setZero();
      for(int s = 0; s < v1.cols(); ++s)
	v1.col(s) = toLogProportionG((vector<Type>)v0.col(s));
      logCM(i) = v1;
    }
    genParSet.logConfusionMatrix = logCM;
  }else{
    // Rcout << "There is NOT a confusion matrix!\n";
    genParSet.logConfusionMatrix = vector<matrix<Type>>(0);
  }
  REPORT(genParSet.logConfusionMatrix);
  PARAMETER_MATRIX(gen_dmScale); genParSet.dmScale = gen_dmScale;
  PARAMETER_MATRIX(gen_muLogP); genParSet.muLogP = gen_muLogP;
  PARAMETER_VECTOR(gen_avgProbPar); genParSet.avgProbPar = gen_avgProbPar;
  PARAMETER_ARRAY(logContamination); // Age x Year x Fleet x (management+Contamination stock); use map to fix things not needed; RE
  PARAMETER_MATRIX(meanContamination);
  PARAMETER_VECTOR(transfRhoContamination);
  PARAMETER_VECTOR(logSdContamination);
 
  PARAMETER_ARRAY(logGst);
  PARAMETER_MATRIX(logGtrip);

  PARAMETER_CMOE_VECTOR(logitArea);


  
  matrix<Type> Parea(stockAreas.rows(), stockAreas.cols());
  Parea.setZero();
  for(int s = 0; s < Parea.cols(); ++s){
    int indx = 0;
    Type ps = 0.0;
    vector<Type> pparsIn = logitArea.col(s);
    vector<Type> pparsOut(Parea.rows());
    pparsOut.setZero();
    for(int a = 0; a < Parea.rows(); ++a){
      if(stockAreas(a,s) == 1){
	if(indx == pparsIn.size()){
	  pparsOut(a) = 1.0;
	  ps += 1.0;
	  break;
	}else{
	  pparsOut(a) = exp(pparsIn(indx++));
	  ps += pparsOut(a);
	}
      }
    }
    if(pparsOut.size() > 0){
      pparsOut /= ps;
      Parea.col(s) = pparsOut;
    }
  }
  REPORT(Parea);


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

  ofall<Type> ofAll(nStocks);

  ////////////////////////////////////
  //////////// Mohn's rho ///////////
  //////////////////////////////////

  if(mohn == 3){ // Use random effects to add correlation
    Type n = sam.dataSets.size();
    int T = sam.dataSets(0).years.size();
    for(int i = 0; i < mohn_fake_obs.size(); ++i){
      if(!isNA(sam.dataSets(0).logobs(i))){
	ans -= n * dnorm(sam.dataSets(0).logobs(i), mohn_fake_obs(i), Type(0.00001), true);
      }else{
       	ans -= n * dnorm(mohn_fake_obs(i), Type(0.0), Type(1.0 / sqrt(2.0 * M_PI)), true);
      }
      
    }
    for(int s = 0; s < n; ++s){ // loop over peels
      Type nobs_s = sam.dataSets(s).logobs.size();
	for(int k = 0; k < nobs_s; ++k){
	  if(!isNA(sam.dataSets(s).logobs(k)))
	    sam.dataSets(s).logobs(k) = mohn_fake_obs(k);
	}
    }
  }

  if(mohn == 2){ // Quadratic form correction

  }


  
  ////////////////////////////////////
  /////////// Missing Obs ///////////
  //////////////////////////////////
  //if(mohn != 3){
    // patch missing single-stock observations with random effects 
    int idxmis = 0;
    for(int s = 0; s < nStocks; ++s)
      for(int i = 0; i < sam.dataSets(s).nobs; ++i){
	if(isNA(sam.dataSets(s).logobs(i))){
	  sam.dataSets(s).logobs(i) = missing(idxmis++);
	  int f = sam.dataSets(s).aux(i,1) - 1;
	  if(sam.dataSets(s).fleetTypes(f)>=90){//==90){
	    ans -= dnorm((Type)missing(idxmis-1),Type(0.0),Type(1.0 / sqrt(2.0 * M_PI)),true);
	  }
	  if(sam.dataSets(s).fleetTypes(f) == 0 && sharedObs.hasSharedObs && sharedObs.keyFleetStock(f,s) == 0)
	    ans -= dnorm((Type)missing(idxmis-1),Type(0.0),Type(1.0 / sqrt(2.0 * M_PI)),true);
	}
      }
    // patch missing vulnerability keys with parameters
    idxmis = 0;
    for(int i = 0; i < sharedObs.keyFleetStock.rows(); ++i)
      for(int j = 0; j < sharedObs.keyFleetStock.cols(); ++j)
	if(isNA(sharedObs.keyFleetStock(i,j))){
	  sharedObs.keyFleetStock(i,j) = invlogit(shared_logitMissingVulnerability(idxmis++));
	}
  
    // patch missing shared observations with random effects 
    idxmis = 0;
    if(sharedObs.hasSharedObs)
      for(int i = 0; i < sharedObs.logobs.size(); ++i)
	if(isNA(sharedObs.logobs(i))){
	  sharedObs.logobs(i) = shared_missingObs(idxmis++);
	}    
  
    // add wide prior for missing random effects, but _only_ when computing ooa residuals
    if(doResiduals){
      Type huge = 10.0;
      for (int i = 0; i < missing.size(); i++)
	ans -= dnorm(missing(i), Type(0.0), huge, true);
      for (int i = 0; i < shared_missingObs.size(); i++)
	ans -= dnorm(shared_missingObs(i), Type(0.0), huge, true);
    } 

    //}

  ////////////////////////////////////
  /////////// Recruitment ///////////
  //////////////////////////////////
  vector<Recruitment<Type> > recruits(nStocks);
  for(int s = 0; s < nStocks; ++s){
    recruits(s) = makeRecruitmentFunction(sam.dataSets(s), sam.confSets(s), paraSets(s));
  }



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

    ans += nllSW(logSWa, ds, cs, ps, sam.forecastSets(s), &of);
    sam.dataSets(s).stockMeanWeight = ds.stockMeanWeight;
    ans += nllCW(logCWa, ds, cs, ps, sam.forecastSets(s), &of);
    sam.dataSets(s).catchMeanWeight = ds.catchMeanWeight;
    ans += nllMO(logitMOa, ds, cs, ps, sam.forecastSets(s), &of);
    sam.dataSets(s).propMat = ds.propMat;
    ans += nllNM(logNMa, ds, cs, ps, sam.forecastSets(s), &of);      
    sam.dataSets(s).natMor = ds.natMor;


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

  ////////////////////////////////////////
  /////////// Prepare forecast //////////
  //////////////////////////////////////

  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    // Calculate forecast
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    // Resize arrays
    prepareForForecast(sam.forecastSets(s), sam.dataSets(s), sam.confSets(s), paraSets(s), logFa, logNa, recruits(s), &of);
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }

  // Update mortalities ahead of time (but after prepare forecast)
  vector<MortalitySet<Type> > mortalities(nStocks);
  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);     
    array<Type> logFa = getArray(logF, s);
    array<Type> lfs = getArray(logitFseason,s);
    mortalities(s) = MortalitySet<Type>(sam.dataSets(s), sam.confSets(s), paraSets(s), logFa, lfs);
    // reportMort((MortalitySet<Type>)mortalities(s),&of);
    MortalitySet<Type> mort = mortalities(s);
    REPORT_F(mort, (&of));
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }


  
  // // Calculate forecast ahead of time
  // for(int s = 0; s < nStocks; ++s){
  //   // Calculate forecast
  //   array<Type> logNa = getArray(logN, s);
  //   array<Type> logFa = getArray(logF, s);
  //   array<Type> logitFSa = getArray(logitFseason, s);
  //   // Resize arrays
  //   sam.forecastSets(s).calculateForecast(logFa,logNa, logitFSa, sam.dataSets(s), sam.confSets(s), paraSets(s), recruits(s), mortalities(s));
  // }
  
  // Overwrite keyLogFsta for observations
  if(sharedObs.hasSharedObs){
    for(int s = 0; s < nStocks; ++s){
      array<Type> logFa = getArray(logF, s);
      confSet cs = sam.confSets(s);
      for(int f = 0; f < cs.keyLogFsta.dim(0); ++f){
	if(sharedObs.keyFleetStock(f,s) == 0){
	  for(int a = 0; a < cs.keyLogFsta.dim(1); ++a){
	    int fnum = sam.confSets(s).keyLogFsta(f,a);	    
	    if(fnum >= 0){
	      for(int y = 0; y < logFa.dim[1]; ++y)
		ans -= dnorm(logFa(fnum,y),Type(0.0),Type(1.0/(2.0*M_PI)),true);
	      sam.confSets(s).keyLogFsta(f,a) = -1;
	    }
	  }
	}
      }
    }
  }  


  vector<bool> hasPH(Xph.size());
  hasPH.setConstant(false);
  //vector<matrix<Type> > phPred(Xph.size());
  vector<array<Type> > phPred(Xph.size());
  for(int s = 0; s < phPred.size(); ++s){
    if(shared_phbeta.col(s).size() > 0){
      vector<Type> v0 = (vector<Type>)(Xph(s).matrix() * shared_phbeta.col(s));
      // if(v0.size() != (maxAgeAll - minAgeAll + 1)*(maxYearAll - minYearAll + 1))
      // 	Rf_error("Wrong size in proportional hazard");
      // phPred(s) = asMatrix(v0, maxAgeAll - minAgeAll + 1, maxYearAll - minYearAll + 1); // Age x stock
      int nAge = maxAgeAll - minAgeAll + 1;
      // int nYear = maxYearAll - minYearAll + 1;
      // int nFleet = v0.size() / nAge / nYear;
      int nFleet = sam.dataSets(s).catchMeanWeight.dim[2];//v0.size() / nAge / nYear;
      int nYear = v0.size() / nAge / nFleet; // Can be less than
      //phPred(s) = asMatrix(v0, nAge, v0.size() / nAge); // Age x year
      vector<int> dd(3);
      dd << nAge, nYear, nFleet;
      array<Type> tmp(v0,dd);//nAge,nYear,nFleet);
      //phPred(s).setdim(dd);
      phPred(s) = tmp;
      hasPH(s) = true;
    }else{
      array<Type> m0(0,0,0);
      m0.setZero();
      phPred(s) = m0;
    }
  }
  //REPORT(phPred);

  // Initial parameter contribution
  for(int s = 0; s < nStocks; ++s){
    if(initLogF.col(s).size() > 0){
      if(shared_logFscale.col(s).cols() == 0){
	for(int a = 0; a < initLogF.col(s).size(); ++a)
	  ans -= dnorm(logF.col(s)(a,0), initLogF.col(s)(a), Type(0.01), true);
      }else{
	ans -= dnorm(shared_logFscale.col(s)(0,0), initLogF.col(s)(0), Type(0.01), true);
      }
    }else if(doResiduals){
      Type huge = 10.0;
      for(int a = 0; a < logF.col(s).rows(); ++a)
	ans -= dnorm(logF.col(s)(a,0), Type(0.0), huge, true);

    }
  }

  // bool overwriteF = false;
  for(int s = 1; s < nStocks; ++s){ // First F is always an RW process
    oftmp<Type> of(this->do_simulate);
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    array<Type> logitFSa = getArray(logitFseason, s);
    data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    confSet cs0 = sam.confSets(0);
    paraSet<Type> ps = paraSets(s);
    forecastSet<Type> fs = sam.forecastSets(s);
    //if(shared_logFscale.col(s).size() == 0 && !hasPH(s)){ // Not using shared logF selectivity or proportional hazard
    if(shared_F_type == -1){ // Same F for all stocks
      matrix<Type> logFs = logF.col(s);
      matrix<Type> logF0 = logF.col(0);
      // Add fake likelihood contribution for unused random effects
      for(int i = 0; i < logFs.cols(); ++i){ // Loop over time
	if(fs.nYears == 0 || fs.forecastYear(i) == 0){
	  for(int j = 0; j < logFs.rows(); ++j){	  
	    ans -= dnorm(logFs(j,i), Type(0.0), Type(1.0 / sqrt(2.0 * M_PI)), true);
	    logFs(j,i) = logF0(j,i);
	  }
	  if(sharedObs.hasSharedObs)
	    for(int f = 0; f < cs.keyLogFsta.dim(0); ++f)
	      for(int a = 0; a < cs.keyLogFsta.dim(1); ++a)
		if(cs0.keyLogFsta(f,a) > (-1))
		  if(sharedObs.keyFleetStock(f,s) == 0)
		    logFs(cs0.keyLogFsta(f,a),i) = R_NegInf;
	}
      }
      // if(fs.nYears > 0){
      // 	array<Type> lFF1 = toArray(logFs);	
      // 	Type n1 = nllF(ds, cs, ps, fs, lFF1, logNa, keepTmp, &of);
      // 	matrix<Type> lFFtmp = logFs.block(0,0,logFs.rows(),fs.preYears);
      // 	array<Type> lFF = toArray(lFFtmp);
      // 	Type n2 = nllF(ds, cs, ps, fs, lFF, logNa, keepTmp, &of);
      // 	ans += n1 - n2;
      // }
       // Overwrite logF
      logF.col(s) = logFs;
    }else if(shared_F_type == 0){ // Nothing shared
      // SAVE FOR LATER
      // ans += nllF(ds, cs, ps, fs, logFa, keepTmp, &of);
      // ofAll.addToReport(of.report,s);
      // moveADREPORT(&of,this,s);
      // // If simulate -> move grab new logF values and move them to the right place!         
      // logF.col(s) = logFa.matrix();
      // Overwrite unused F with -Inf
      if(sharedObs.hasSharedObs){
	// matrix<Type> logFs = logF.col(s);     
	// for(int i = 0; i < logFs.cols(); ++i){ 
	//   for(int f = 0; f < cs.keyLogFsta.dim(0); ++f)
	//     for(int a = 0; a < cs.keyLogFsta.dim(1); ++a)
	//       //if(cs0.keyLogFsta(f,a) > (-1))
	//       if(sharedObs.keyFleetStock(f,s) == 0 & sharedObs.fleetType(f) == 0){
	// 	  ans -= dnorm(logFs(f,a), Type(0.0), Type(1.0 / sqrt(2.0 * M_PI)), true);
	// 	  logFs(cs0.keyLogFsta(f,a),i) = R_NegInf;		  
	//       }
	// }
	// logF.col(s) = logFs;
      }
    }else if(shared_F_type == 1){ // Scaling by vector AR(1)
      matrix<Type> logF0 = logF.col(0);
      matrix<Type> logFs = logF.col(s);
      for(int i = 0; i < logFs.cols(); ++i){ // Loop over time
	if(fs.nYears == 0 || fs.forecastYear(i) == 0){
	  for(int a = 0; a < logF0.rows(); ++a){ // Loop over F ages
	    Type rho = toInterval(shared_lfsRho(s-1),Type(0.0),Type(1.0),Type(2.0));
	    Type mu = shared_lfsMean(a,s-1);
	    if(i > 0)
	      mu += rho * (shared_logFscale.col(s)(a,i-1) - shared_lfsMean(a,s-1));
	    Type sd = exp(shared_lfsSd(s-1));
	    if(i == 0)
	      sd /= sqrt(1.0 - rho * rho);
	    ans -= dnorm(shared_logFscale.col(s)(a,i), mu, sd, true);
	    ans -= dnorm(logFs(a,i), Type(0.0), Type(1.0 / sqrt(2.0 * M_PI)), true);
	    logFs(a,i) = logF0(a,i) + shared_logFscale.col(s)(a,i);
	  }
	  if(sharedObs.hasSharedObs)
	    for(int f = 0; f < cs.keyLogFsta.dim(0); ++f)
	      for(int a = 0; a < cs.keyLogFsta.dim(1); ++a)
		if(cs0.keyLogFsta(f,a) > (-1))
		  if(sharedObs.keyFleetStock(f,s) == 0)
		    logFs(cs0.keyLogFsta(f,a),i) = R_NegInf;
	}else if(fs.nYears > 0 && fs.forecastYear(i) > 0){
	  array<Type> lFF1 = toArray(logFs);
	  fs.updateForecast(CppAD::Integer(fs.forecastYear(i)) - 1, lFF1, logNa, logitFSa, ds, cs, ps, recruits(s), mortalities(s), this->do_simulate);
	}
      }
      // if(fs.nYears > 0){
      // 	array<Type> lFF1 = toArray(logFs);	
      // 	Type n1 = nllF(ds, cs, ps, fs, lFF1, logNa, keepTmp, &of);
      // 	matrix<Type> lFFtmp = logFs.block(0,0,logFs.rows(),fs.preYears);
      // 	array<Type> lFF = toArray(lFFtmp);
      // 	Type n2 = nllF(ds, cs, ps, fs, lFF, logNa, keepTmp, &of);
      // 	ans += n1 - n2;
      // }
      logF.col(s) = logFs;
    }else{			// Any combination of scaling by parametric function and scalar RW
      matrix<Type> logF0 = logF.col(0);
      matrix<Type> logFs = logF.col(s);      
      for(int i = 0; i < logFs.cols(); ++i){ // Loop over time
	if(fs.nYears == 0 || fs.forecastYear(i) == 0){
	  Type slfs = 0.0;
	  if(shared_F_type != 3){	// If not pure parametric
	    if(shared_F_type == 2 || shared_F_type == 4){ // AR in time
	      Type rho = toInterval(shared_lfsRho(s-1),Type(0.0),Type(1.0),Type(2.0));
	      Type mu = shared_lfsMean(0,s-1);
	      if(i > 0)
		mu += rho * (shared_logFscale.col(s)(0,i-1) - shared_lfsMean(0,s-1));	  
	      Type sd = exp(shared_lfsSd(s-1));
	      if(i == 0)
		sd /= sqrt(1.0 - rho * rho);
	      ans -= dnorm(shared_logFscale.col(s)(0,i), mu, sd, true);
	    }else{ // 5 or 6: RW in time
	      Type sd = exp(shared_lfsSd(s-1));
	      if(i == 0){
		ans -= dnorm(shared_logFscale.col(s)(0,i), shared_lfsMean(0,s-1), Type(0.01), true);
	      }else{
		ans -= dnorm(shared_logFscale.col(s)(0,i), shared_logFscale.col(s)(0,i-1), sd, true);
	      }
	    }
	    slfs = shared_logFscale.col(s)(0,i);
	  }
	  vector<Type> lfPred = logF0.col(i);
	  vector<Type> lfAdd(lfPred.size()); lfAdd.setZero();
	  vector<Type> lfNum(lfPred.size()); lfNum.setZero();
	  vector<int> unused(lfPred.size()); unused.setZero();
	  if(shared_F_type != 2 && shared_F_type != 5  && hasPH(s)){
	    for(int f = 0; f < cs.keyLogFsta.dim(0); ++f){
	      for(int a = 0; a < cs.keyLogFsta.dim(1); ++a){
		if(cs0.keyLogFsta(f,a) > (-1)){
		  lfAdd(cs0.keyLogFsta(f,a)) += (phPred(s).col(f))(a,std::min(i,(int)phPred(s).col(f).cols()-1));
		  lfNum(cs0.keyLogFsta(f,a)) += 1.0;
		  if(sharedObs.hasSharedObs)
		    if(sharedObs.keyFleetStock(f,s) == 0){
		      unused(cs0.keyLogFsta(f,a)) = 1;
		    }
		}
	      }
	    }
	    for(int j = 0; j < lfPred.size(); ++j){
	      if(lfNum(j) > 0)
		lfPred(j) += lfAdd(j) / lfNum(j);
	      ans -= dnorm(logFs(j,i), Type(0.0), Type(1.0 / sqrt(2.0 * M_PI)), true);
	      if(unused(j)){
		logFs(j,i) = R_NegInf;
	      }else{
		logFs(j,i) = lfPred(j) + slfs;
	      }
	    }
	  }
	  // Implement simulation of logF with shared selectivity(?)
	  array<Type> lFF1 = toArray(logFs);
	  mortalities(s).updateYear(ds, cs, ps, lFF1, logitFSa,i);
	}else if(fs.nYears > 0 && fs.forecastYear(i) > 0){
	  // fake contribution from unused shared_logFscale
	  // for(int j = 0; j < logFs.rows(); ++j){
	  //   Rcout << "(" << j << "," << i << "): " << logFs(j,i) << ", " << logFa(j,i) << "\n";
	  //   logFs(j,i) = logFa(j,i);
	  // }
	  array<Type> lFF1 = toArray(logFs);
	  fs.updateForecast(CppAD::Integer(fs.forecastYear(i)) - 1, lFF1, logNa, logitFSa, ds, cs, ps, recruits(s), mortalities(s), this->do_simulate);
	  mortalities(s).updateYear(ds, cs, ps, lFF1, logitFSa,i);
	  ans -= dnorm(shared_logFscale.col(s)(0,i), Type(0.0), Type(1.0 / sqrt(2.0 * M_PI)), true);
	  
	}
      }      
      // if(fs.nYears > 0){
      // 	array<Type> lFF1 = toArray(logFs);	
      // 	Type n1 = nllF(ds, cs, ps, fs, lFF1, logNa, keepTmp, &of);
      // 	matrix<Type> lFFtmp = logFs.block(0,0,logFs.rows(),fs.preYears);
      // 	array<Type> lFF = toArray(lFFtmp);
      // 	Type n2 = nllF(ds, cs, ps, fs, lFF, logNa, keepTmp, &of);
      // 	ans += n1 - n2;
      // }
      logF.col(s) = logFs;
      //overwriteF = true;
    }
  }

  // if(overwriteF){
  for(int s = 0; s < logF.cols(); ++s){
    oftmp<Type> of(this->do_simulate);
    matrix<Type> logFs = logF.col(s);
    REPORT_F(logFs,(&of));
    ADREPORT_F(logFs,(&of));
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }
  // }

  if(shared_Fseason_type == 1){
    for(int s = 0; s < nStocks; ++s){      
      array<Type> logitFSa = getArray(logitFseason, s);
      logitFseason.col(s) = logitFseason.col(0);
      for(int i = 0; i < logitFseason.col(0).size(); ++i)
	ans -= dnorm(logitFSa[i], Type(0.0), Type(1.0 / sqrt(2.0 * M_PI)), true);
    }
  }


   ////////////////////////////////
  ////////// Season /////////////
  //////////////////////////////

   for(int s = 0; s < nStocks; ++s){
     oftmp<Type> of(this->do_simulate);
     array<Type> logitFSa = getArray(logitFseason, s);
     data_indicator<vector<Type>,Type> keepTmp = keep(s);
     dataSet<Type> ds = sam.dataSets(s);
     confSet cs = sam.confSets(s);
     paraSet<Type> ps = paraSets(s);
     if(shared_Fseason_type == 0 || s == 0){
       ans += nllSeason(ds, cs, ps, sam.forecastSets(s), logitFSa, keepTmp, &of);
       ofAll.addToReport(of.report,s);
       moveADREPORT(&of,this,s);
     // If simulate -> move grab new logF values and move them to the right place!         
       logitFseason.col(s) = logitFSa;
     }
   }


   ////////////////////////////////////////
   /////////// Mortality /////////////////
   //////////////////////////////////////
   //vector<MortalitySet<Type> > mortalities(nStocks);
   // Update mortalities
   for(int s = 0; s < nStocks; ++s){
     oftmp<Type> of(this->do_simulate);     
     array<Type> logFa = getArray(logF, s);
     array<Type> lfs = getArray(logitFseason,s);
     mortalities(s) = MortalitySet<Type>(sam.dataSets(s), sam.confSets(s), paraSets(s), logFa, lfs);
     // reportMort((MortalitySet<Type>)mortalities(s),&of);
     MortalitySet<Type> mort = mortalities(s);
     REPORT_F(mort, (&of));
     ofAll.addToReport(of.report,s);
     moveADREPORT(&of,this,s);
   }



  ////////////////////////////////////////
  /////////// Calculate forecast ////////
  //////////////////////////////////////

   // Update calculated forecast
   for(int s = 0; s < nStocks; ++s){
     // Calculate forecast
     array<Type> logNa = getArray(logN, s);
     array<Type> logFa = getArray(logF, s);
     array<Type> logitFSa = getArray(logitFseason, s);
     // Resize arrays
     sam.forecastSets(s).calculateForecast(logFa,logNa, logitFSa, sam.dataSets(s), sam.confSets(s), paraSets(s), recruits(s), mortalities(s));
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
    //if(shared_logFscale.col(s).size() == 0 && !hasPH(s)){ // Not using shared logF selectivity or proportional hazard
    if(shared_F_type == 0 || s == 0){
      ans += nllF(ds, cs, ps, sam.forecastSets(s), logFa, logNa, keepTmp, &of);
      ofAll.addToReport(of.report,s);
      moveADREPORT(&of,this,s);
      // If simulate -> move grab new logF values and move them to the right place!
      // Does this ruin things???
      //logF.col(s) = asMatrix<Type>((vector<Type>)logFa.vec(),logFa.dim[0],logFa.dim[1]);
      logF.col(s) = logFa.matrix();
    }else if(shared_F_type != 0 && s > 0){
      // Get likelihood contribution of forecast
      // if(sam.forecastSets(s).nYears > 0){
      // 	Rcout << s << ": " << logFa.dim << "\n";
      // 	Rcout << "ny: " << sam.forecastSets(s).nYears << "\n";
      // 	Rcout << "py: " << sam.forecastSets(s).preYears << "\n";
      // 	array<Type> lFF1 = logFa;	
      // 	Type n1 = nllF(ds, cs, ps, sam.forecastSets(s), lFF1, logNa, keepTmp, &of);
      // 	matrix<Type> lFFtmp = logFa.matrix().block(0,0,logFa.dim[0],sam.forecastSets(s).preYears);
      // 	array<Type> lFF = toArray(lFFtmp);
      // 	Type n2 = nllF(ds, cs, ps, sam.forecastSets(s), lFF, logNa, keepTmp, &of);
      // 	ans += n1 - n2;
      // }
      // Update other stocks if simulating historical values
      SIMULATE{
	GetRNGstate();
	if(cs.simFlag(0)==0){
	  if(shared_F_type == 1){ // Scaling by vector AR(1)
	    matrix<Type> logF0 = logF.col(0);
	    matrix<Type> logFs = logF.col(s);
	    for(int i = 0; i < logFs.cols(); ++i){ // Loop over time
	      if(sam.forecastSets(s).nYears == 0 || sam.forecastSets(s).forecastYear(i) == 0){
		for(int a = 0; a < logF0.rows(); ++a){ // Loop over F ages
		  Type rho = toInterval(shared_lfsRho(s-1),Type(0.0),Type(1.0),Type(2.0));
		  Type mu = shared_lfsMean(a,s-1);
		  if(i > 0)
		    mu += rho * (shared_logFscale.col(s)(a,i-1) - shared_lfsMean(a,s-1));
		  Type sd = exp(shared_lfsSd(s-1));
		  if(i == 0)
		    sd /= sqrt(1.0 - rho * rho);		  
		  shared_logFscale.col(s)(a,i) = rnorm(1,mu, sd)[0];
		  logFs(a,i) = logF0(a,i) + shared_logFscale.col(s)(a,i);
		}
	      }
	    }
	    logF.col(s) = logFs;
	    REPORT_F(logFs,(&of));
	    ADREPORT_F(logFs,(&of));
	    ofAll.addToReport(of.report,s);
	    moveADREPORT(&of,this,s);
	  }else{			// Any combination of scaling by parametric function and scalar RW
	    matrix<Type> logF0 = logF.col(0);
	    matrix<Type> logFs = logF.col(s);
	    confSet cs0 = sam.confSets(0);
	    for(int i = 0; i < logFs.cols(); ++i){ // Loop over time
	      if(sam.forecastSets(s).nYears == 0 || sam.forecastSets(s).forecastYear(i) == 0){
		Type slfs = 0.0;
		if(shared_F_type != 3){	// If not pure parametric
		  if(shared_F_type == 2 || shared_F_type == 4){ // AR in time
		    Type rho = toInterval(shared_lfsRho(s-1),Type(0.0),Type(1.0),Type(2.0));
		    Type mu = shared_lfsMean(0,s-1);
		    if(i > 0)
		      mu += rho * (shared_logFscale.col(s)(0,i-1) - shared_lfsMean(0,s-1));	  
		    Type sd = exp(shared_lfsSd(s-1));
		    if(i == 0)
		      sd /= sqrt(1.0 - rho * rho);
		    shared_logFscale.col(s)(0,i) = rnorm(1,mu, sd)[0];
		  }else{ // 5 or 6: RW in time
		    Type sd = exp(shared_lfsSd(s-1));
		    if(i == 0){
		      shared_logFscale.col(s)(0,i) = rnorm(1,shared_lfsMean(0,s-1), Type(0.01))[0];
		    }else{
		      shared_logFscale.col(s)(0,i) = rnorm(1,shared_logFscale.col(s)(0,i-1), sd)[0];
		    }
		  }
		  slfs = shared_logFscale.col(s)(0,i);
		}
		vector<Type> lfPred = logF0.col(i);
		vector<Type> lfAdd(lfPred.size()); lfAdd.setZero();
		vector<Type> lfNum(lfPred.size()); lfNum.setZero();
		vector<int> unused(lfPred.size()); unused.setZero();
		if(shared_F_type != 2 && shared_F_type != 5  && hasPH(s)){
		  for(int f = 0; f < cs.keyLogFsta.dim(0); ++f)
		    for(int a = 0; a < cs.keyLogFsta.dim(1); ++a)
		      if(cs0.keyLogFsta(f,a) > (-1)){
			lfAdd(cs0.keyLogFsta(f,a)) += phPred(s)(a,std::min(i,(int)phPred(s).cols()-1));
			lfNum(cs0.keyLogFsta(f,a)) += 1.0;
			if(sharedObs.hasSharedObs)
			  if(sharedObs.keyFleetStock(f,s) == 0){
			    unused(cs0.keyLogFsta(f,a)) = 1;
			  }
		      }
		  for(int j = 0; j < lfPred.size(); ++j){
		    if(lfNum(j) > 0)
		      lfPred(j) += lfAdd(j) / lfNum(j);
		    //ans -= dnorm(logFs(j,i), Type(0.0), Type(1.0 / sqrt(2.0 * M_PI)), true);
		    if(unused(j)){
		      logFs(j,i) = R_NegInf;
		    }else{
		      logFs(j,i) = lfPred(j) + slfs;
		    }
		  }
		}
		// Implement simulation of logF with shared selectivity(?)	 
	      }
	    }
	    logF.col(s) = logFs;
	    REPORT_F(logFs,(&of));
	    ADREPORT_F(logFs,(&of));
	    ofAll.addToReport(of.report,s);
	    moveADREPORT(&of,this,s);
	  }
	}
	PutRNGstate();
      }
    }
   }
   // Update mortalities if simulating historical values
   SIMULATE{
     for(int s = 0; s < logF.cols(); ++s){
       oftmp<Type> of(this->do_simulate);
       matrix<Type> logFs = logF.col(s);
       REPORT_F(logFs,(&of));
       ADREPORT_F(logFs,(&of));
       ofAll.addToReport(of.report,s);
       moveADREPORT(&of,this,s);
     }
     for(int s = 0; s < nStocks; ++s){
       oftmp<Type> of(this->do_simulate);     
       array<Type> logFa = getArray(logF, s);
       array<Type> lfs = getArray(logitFseason,s);
       mortalities(s) = MortalitySet<Type>(sam.dataSets(s), sam.confSets(s), paraSets(s), logFa, lfs);
       // reportMort((MortalitySet<Type>)mortalities(s),&of);
       MortalitySet<Type> mort = mortalities(s);
       REPORT_F(mort, (&of));
       ofAll.addToReport(of.report,s);
       moveADREPORT(&of,this,s);
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
	for(int f = 0; f < conf.keyLogFsta.dim(0); ++f)
	  if(conf.keyLogFsta(f,a-1)>(-1))
	    m -= exp(logFa(conf.keyLogFsta(f,a-1),0));	
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
  REPORT(ncov);

    ////////////////////////////////////
    ////// CALCULATE CONTRIBUTION //////
    ////////////////////////////////////


    // add wide prior for first state, but _only_ when computing osa residuals
    if(doResiduals){
      Type huge = 10.0;
      for(int s = 0; s < nAreas; ++s){
	if(initLogN.col(s).size() == 0)
	  for (int i = 0; i < logN.col(s).rows(); i++){
	    ans -= dnorm((logN.col(s))(i, 0), Type(0.0), huge, true);
	  }
      }
    } 

 

   density::MVNORM_t<Type> neg_log_densityN(ncov);


   vector<MVMIX_t<Type> > neg_log_densityF(nAreas);
   for(int s = 0; s < nAreas; ++s){
     array<Type> logFa = getArray(logF, s);
     dataSet<Type> ds = sam.dataSets(s);
     confSet cs = sam.confSets(s);
     paraSet<Type> ps = paraSets(s);     
     matrix<Type> fvar = get_fvar(ds, cs, ps, logFa);
     if(sam.forecastSets(s).FEstCov.cols() > 0)
       fvar = sam.forecastSets(s).FEstCov;
     neg_log_densityF(s) = MVMIX_t<Type>(fvar,Type(cs.fracMixF));
   }

   // Loop over time
   for(int yall = 0; yall < maxYearAll - minYearAll + 1; ++yall){
     // Handle simulation of F for forecast and HCR
     SIMULATE{
       GetRNGstate();
       // Simulate new F for forecast and HCR
       for(int s = 0; s < nAreas; ++s){
	 array<Type> logNa = getArray(logN, s);
	 array<Type> logFa = getArray(logF, s);
	 array<Type> logitFSa = getArray(logitFseason, s);
	 //forecastSet<Type> fc = sam.forecastSets(s); // While convenient with shorthand, this did not keep updates made by updateForecast!
	 if(sam.forecastSets(s).nYears > 0){
	   dataSet<Type> ds = sam.dataSets(s);
	   confSet cs = sam.confSets(s);
	   paraSet<Type> ps = paraSets(s);
	   int y = yall - CppAD::Integer(ds.years(0) - minYearAll);
	   int fcOffset = 0;
	   if(sam.forecastSets(s).preYears > 0)
	     fcOffset = ds.noYears - sam.forecastSets(s).preYears;
	   if(y > 0 && y < ds.noYears + sam.forecastSets(s).nYears - fcOffset){
	     // matrix<Type> fvar = get_fvar(ds, cs, ps, logFa);
	     // if(sam.forecastSets(s).FEstCov.cols() > 0)
	     //   fvar = sam.forecastSets(s).FEstCov;
	     // MVMIX_t<Type> neg_log_densityF(fvar,Type(cs.fracMixF));
	     //int nYears = sam.forecastSets(s).nYears;
	     // int fi = y - sam.forecastSets(s).preYears;
	     int fi = CppAD::Integer(sam.forecastSets(s).forecastYear(y)) - 1;
	     //int fi = y - sam.forecastSets(s).forecastYear.size() + nYears;
	     // Update forecast
	     if(fi >= 0){
	       // SET sam.forecastSets(s).nvar HERE!!! SAM CODE WILL NOT WORK WITH CORRELATION IN N
	       sam.forecastSets(s).updateForecast(fi, logFa, logNa, logitFSa, ds, cs, ps, recruits(s), mortalities(s), this->do_simulate);
	       // Simulate F
	       // int forecastIndex = CppAD::Integer(dat.forecast.forecastYear(i))-1;
	       if(sam.forecastSets(s).simFlag(0) == 0){
		 Type timeScale = exp(sam.forecastSets(s).forecastCalculatedLogSdCorrection(fi));
		 if(sam.forecastSets(s).fsdTimeScaleModel(fi) == sam.forecastSets(s).fixedDeviation){
		   logFa.col(y) = (vector<Type>)sam.forecastSets(s).forecastCalculatedMedian.col(fi);
		 }else{
		     logFa.col(y) = (vector<Type>)sam.forecastSets(s).forecastCalculatedMedian.col(fi) + neg_log_densityF(s).simulate() * timeScale;
		 }
		 logF.col(s).col(y) = logFa.col(y);
		 mortalities(s).updateYear(ds, cs, ps, logFa, logitFSa,y);
	       }	   
	     }	 
	   }
	 }
       }
       PutRNGstate();
     }else{
       // Get likelihood contribution of forecast for stock 2++
       for(int s = 1; s < nAreas; ++s){
	 array<Type> logNa = getArray(logN, s);
	 array<Type> logFa = getArray(logF, s);
	 array<Type> logitFSa = getArray(logitFseason, s);
	 if(sam.forecastSets(s).nYears > 0){
	   dataSet<Type> ds = sam.dataSets(s);
	   confSet cs = sam.confSets(s);
	   paraSet<Type> ps = paraSets(s);
	   int y = yall - CppAD::Integer(ds.years(0) - minYearAll);
	   int fcOffset = 0;
	   if(sam.forecastSets(s).preYears > 0)
	     fcOffset = ds.noYears - sam.forecastSets(s).preYears;
	   if(y > 0 && y < ds.noYears + sam.forecastSets(s).nYears - fcOffset){
	     // matrix<Type> fvar = get_fvar(ds, cs, ps, logFa);
	     // if(sam.forecastSets(s).FEstCov.cols() > 0)
	     //   fvar = sam.forecastSets(s).FEstCov;
	     // MVMIX_t<Type> neg_log_densityF(fvar,Type(cs.fracMixF));
	     //int nYears = sam.forecastSets(s).nYears;
	     // int fi = y - sam.forecastSets(s).preYears;
	     int fi = CppAD::Integer(sam.forecastSets(s).forecastYear(y)) - 1;
	     //int fi = y - sam.forecastSets(s).forecastYear.size() + nYears;
	     // Update forecast
	     if(fi >= 0){
	       // SET sam.forecastSets(s).nvar HERE!!! SAM CODE WILL NOT WORK WITH CORRELATION IN N
	       sam.forecastSets(s).updateForecast(fi, logFa, logNa, logitFSa, ds, cs, ps, recruits(s), mortalities(s), this->do_simulate); // Should not be necessary??
	       logF.col(s).col(y) = logFa.col(y);
	       mortalities(s).updateYear(ds, cs, ps, logFa, logitFSa,y);	       
	       //  F contribution
	       // int forecastIndex = CppAD::Integer(dat.forecast.forecastYear(i))-1;
	       Type timeScale = exp(sam.forecastSets(s).forecastCalculatedLogSdCorrection(fi));
	       ans += neg_log_densityF(s)((logFa.col(y) - (vector<Type>)sam.forecastSets(s).forecastCalculatedMedian.col(fi)) / timeScale) + log(timeScale) * Type(logFa.dim[0]);
	     }	 
	   }
	 }
	 logF.col(s) = logFa.matrix();
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
    vector<Type> logVarScale(predN.size());
    logVarScale.setZero();
    // Loop over stocks
    for(int s = 0; s < nAreas; ++s){
      dataSet<Type> ds = sam.dataSets(s);
      confSet cs = sam.confSets(s);
      paraSet<Type> ps = paraSets(s);
      int ageOffset = sam.confSets(s).minAge - minAgeAll;
      int y = yall - CppAD::Integer(ds.years(0) - minYearAll);
      
      array<Type> logNa = getArray(logN, s);
      array<Type> logFa = getArray(logF, s);
      int fcOffset = 0;
	if(sam.forecastSets(s).preYears > 0)
	  fcOffset = ds.noYears - sam.forecastSets(s).preYears;
      if(y > 0 && y < ds.noYears + sam.forecastSets(s).nYears - fcOffset){
	vector<Type> predNnz = predNFun(ds, cs, ps, logNa, logFa, recruits(s), mortalities(s), (int)y);
	keepN.segment(s * nages + ageOffset,predNnz.size()) = 1.0;
	predN.segment(s * nages + ageOffset,predNnz.size()) = predNnz;
	newN.segment(s * nages + ageOffset,predNnz.size()) = logNa.col(y);
	lastN.segment(s * nages + ageOffset,predNnz.size()) = logNa.col(y-1);	
      }
    }
    predN -= AlphaCon * lastN;
    vector<Type> Nscale(predN.size());
    Nscale.setConstant(1.0);
   if(keepN.sum()>0){
      // forecast correction to recruitment
       // Nscale += 1.0;
      // vector<Type> predNTmp = predN;

      for(int s = 0; s < nAreas; ++s){
	dataSet<Type> ds = sam.dataSets(s);
	confSet cs = sam.confSets(s);
	paraSet<Type> ps = paraSets(s);
	array<Type> logNa = getArray(logN, s);
	array<Type> logFa = getArray(logF, s);
	int ageOffset = sam.confSets(s).minAge - minAgeAll;
	int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	int fcOffset = 0;
	if(sam.forecastSets(s).preYears > 0)
	  fcOffset = ds.noYears - sam.forecastSets(s).preYears;
	if(y > 0 &&
	   y < ds.noYears + sam.forecastSets(s).nYears - fcOffset &&
	   sam.forecastSets(s).nYears > 0 &&
	   sam.forecastSets(s).forecastYear(y) > 0 &&
	   sam.forecastSets(s).recModel(CppAD::Integer(sam.forecastSets(s).forecastYear(y))-1) == sam.forecastSets(s).useIID &&
	   sam.forecastSets(s).forecastYear(y) > 0){
	  Nscale(s * nages + ageOffset) = sqrt(sam.forecastSets(s).logRecruitmentVar(CppAD::Integer(sam.forecastSets(s).forecastYear(y))-1)) / sqrt(ncov(s * nages + ageOffset,s * nages + ageOffset));
	  predN(s * nages + ageOffset) = sam.forecastSets(s).logRecruitmentMedian(CppAD::Integer(sam.forecastSets(s).forecastYear(y))-1);
	}
	if(ps.recVarScalePar.size() > 0){
	  Type sigmaR = sqrt(ncov(s * nages + ageOffset,s * nages + ageOffset));
	  // Type log_kappa2 = sigmaR * sigmaR;
	  // First is for recruitment
	  Type logThisSSB=Type(R_NegInf);
	  if((y-cs.minAge)>=0){
	    // logThisSSB=ssbi(dat,conf,logN,logF,mort,i-conf.minAge, true);    
	    logThisSSB=erbi(ds,cs,ps,logNa,logFa,mortalities(s),y-cs.minAge, true);    
	  }else{
	    // logThisSSB=ssbi(dat,conf,logN,logF,mort,0, true); // use first in beginning       
	    logThisSSB=erbi(ds,cs,ps,logNa,logFa,mortalities(s),0, true); // use first in beginning       
	  }
	  Type xx = 0.0;
	  Type xv = 1.0;
	  for(int qq = 0; qq < ps.recVarScalePar.size(); ++qq){
	    xv *=  (predN(s * nages + ageOffset) - logThisSSB);
	    xx += (exp( ps.recVarScalePar(qq) )-1) * xv;
	    xx += (( ps.recVarScalePar(qq) )) * xv;
	  }
	  // xx += ps.recVarScalePar(0) * logThisSSB;
	  // if(ps.recVarScalePar.size() > 1)
	  //   xx += ps.recVarScalePar(1) * predN(s * nages + ageOffset);
	  
	  Type k = log(sigmaR) + xx;
	  logVarScale(s * nages + ageOffset) = log(sqrt(log(exp(k)+Type(1))) / sigmaR);	  	 
	}
      }
      Nscale *= exp(logVarScale);
      ans+=neg_log_densityN((newN-predN) / Nscale, keepN) + (keepN * log(Nscale)).sum();	  
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
      GetRNGstate();
      // 1) Check if any simFlags are 0
      matrix<Type> NscaleSim(Nscale.size(),Nscale.size());
      NscaleSim.setZero();
      NscaleSim.diagonal() = Nscale;
      matrix<Type> ncovSim = NscaleSim * ncov * NscaleSim; // No need to transpose a diagonal matrix
      bool doSim = false;
      bool isForecast = false;
      // Check if at least one stock should be simulated
      for(int s = 0; s < nAreas; ++s){
	int nYears = sam.forecastSets(s).nYears;
	if(nYears > 0){
	  int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	  int fcOffset = 0;
	   if(sam.forecastSets(s).preYears > 0)
	     fcOffset = sam.dataSets(s).noYears - sam.forecastSets(s).preYears;
	  //int fi = y - sam.forecastSets(s).preYears;// y - sam.forecastSets(s).forecastYear.size() + nYears;
	  int fi = -100;	  
	  if(y > 0 &&
	     y < sam.forecastSets(s).forecastYear.size() && //sam.dataSets(s).noYears + sam.forecastSets(s).nYears - fcOffset &&
	   sam.forecastSets(s).nYears > 0 &&
	   sam.forecastSets(s).forecastYear(y) > 0 &&
	     //sam.forecastSets(s).recModel(CppAD::Integer(sam.forecastSets(s).forecastYear(y))-1) != sam.forecastSets(s).asRecModel &&
	     sam.forecastSets(s).forecastYear(y) > 0){
	    isForecast = TRUE;
	    fi = CppAD::Integer(sam.forecastSets(s).forecastYear(y)) - 1;
	  }
	  if(fi >= 0 && sam.forecastSets(s).simFlag(1) == 0){
	    if(fi == 0 && sam.forecastSets(s).useModelLastN){
	      //doSim = false;
	    }else{
	      doSim = true;
	    }
	    break;
	  }else{
	    //doSim = false;
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
	vector<int> notCondOn(ncovSim.cols());
	notCondOn.setZero();
	for(int i = 0; i < notCondOn.size(); ++i){
	  int s = (int)i / (int)nages;
	  int a = (int)i % (int)nages;
	  int nYears = sam.forecastSets(s).nYears;
	  if(nYears > 0){
	    int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	    //int fi = y - sam.forecastSets(s).preYears; // y - sam.forecastSets(s).forecastYear.size() + nYears;
	    //int fi = CppAD::Integer(sam.forecastSets(s).forecastYear(y)) - 1;
	    int fcOffset = 0;
	    if(sam.forecastSets(s).preYears > 0)
	      fcOffset = sam.dataSets(s).noYears - sam.forecastSets(s).preYears;
	    int fi = -100;
	    if(y > 0 &&
	       y <  sam.forecastSets(s).forecastYear.size() && //sam.dataSets(s).noYears + sam.forecastSets(s).nYears - fcOffset &&
	       sam.forecastSets(s).nYears > 0 &&
	       sam.forecastSets(s).forecastYear(y) > 0 &&
	       //sam.forecastSets(s).recModel(CppAD::Integer(sam.forecastSets(s).forecastYear(y))-1) != sam.forecastSets(s).asRecModel &&
	       sam.forecastSets(s).forecastYear(y) > 0){
	      fi = CppAD::Integer(sam.forecastSets(s).forecastYear(y)) - 1;
	    }
	      if((sam.confSets(s).simFlag(1) == 0 ||
		  (fi >= 0 && sam.forecastSets(s).simFlag(1) == 0)) &&
		 keepN(i) == 1){
		if(!(fi == 0 && sam.forecastSets(s).useModelLastN))
		  notCondOn(i) = 1;
	      }
	  }else if(sam.confSets(s).simFlag(1) == 0 && keepN(i) == 1){
	    notCondOn(i) = 1;
	  }	  
	}
	int nAll = ncovSim.cols();
	int nNotCond = notCondOn.sum();
	int nCond = nAll - nNotCond;
	if(nCond == 0){ // Simulate all stocks
	  vector<Type> noiseN = density::MVNORM(ncovSim).simulate();
	  vector<Type> simRes = predN + noiseN;
	  REPORT(noiseN);
	  // Insert into logN
	  int k1 = 0;
	  for(int i = 0; i < simRes.size(); ++i){	    
	    if(notCondOn(i) == 1){
	      int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
	      // Age index = age - minAge
	      int a = i % nages;
	      int ageOffset = sam.confSets(s).minAge - minAgeAll;
	      int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	      Type logRec = logN.col(s)(0,y);
	      logN.col(s)(a - ageOffset,y) = simRes(k1++);
	      if(a - ageOffset == 0 && sam.confSets(s).simKeepRec){// && isForecast)
		logN.col(s)(a - ageOffset,y) = logRec;
	      }
	    }
	  }
	  // Handle age zero recruitment
	    bool hasZeroRec = false;
	  k1 = 0; int k2 = 0;
	  for(int i = 0; i < notCondOn.size(); ++i){
	    int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
	    // Age index = age - minAge
	    int a = i % nages;
	    int ageOffset = sam.confSets(s).minAge - minAgeAll;
	    int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	    if(sam.confSets(s).minAge == 0 &&
	       sam.forecastSets(s).recModel(CppAD::Integer(sam.forecastSets(s).forecastYear(y))-1) == sam.forecastSets(s).asRecModel ){
	      hasZeroRec = true;
	      array<Type> logNa = getArray(logN, s);
	      array<Type> logFa = getArray(logF, s);
	      dataSet<Type> ds = sam.dataSets(s);
	      confSet cs = sam.confSets(s);
	      paraSet<Type> ps = paraSets(s);
	      // Overwrite (only) recruitment
	      if(a - ageOffset == 0)
		predN(i) = predNFun(ds, cs, ps, logNa, logFa, recruits(s), mortalities(s), (int)y)(0);	     
	    }
	  }
	  // Update simres and insert
	  if(hasZeroRec){
	    simRes = predN + noiseN;
	    k1 = 0;
	    for(int i = 0; i < notCondOn.size(); ++i){
	      if(notCondOn(i) == 1){
		int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
		// Age index = age - minAge
		int a = i % nages;
		int ageOffset = sam.confSets(s).minAge - minAgeAll;
		int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
		Type logRec = logN.col(s)(0,y);
		logN.col(s)(a - ageOffset,y) = simRes(k1++);
		if(a - ageOffset == 0 && sam.confSets(s).simKeepRec && isForecast){
		  logN.col(s)(a - ageOffset,y) = logRec;
		}
	      }
	    }
	  }
	}else if(nNotCond > 0){ // Only do this if there are any not conditional
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
	      Sigma_NN(i,j) = ncovSim(ccond(i),ccond(j));

	  for(int i = 0; i < cond.size(); ++i)
	    for(int j = 0; j < cond.size(); ++j)
	      Sigma_CC(i,j) = ncovSim(cond(i),cond(j));
	  matrix<Type> Sigma_CC_inv = Sigma_CC.inverse();

	  for(int i = 0; i < ccond.size(); ++i)
	    for(int j = 0; j < cond.size(); ++j)
	      Sigma_NC(i,j) = ncovSim(ccond(i),cond(j));
	  
	  matrix<Type> meanCorrectionMat = Sigma_NC * Sigma_CC_inv;

	  for(int i = 0; i < cond.size(); ++i)
	    for(int j = 0; j < ccond.size(); ++j)
	      Sigma_CN(i,j) = ncovSim(cond(i),ccond(j));

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
	  vector<Type> noiseN = density::MVNORM(newSigma).simulate();
	  vector<Type> simRes = muNew + (vector<Type>)(meanCorrectionMat * avec) + noiseN;
	  REPORT(noiseN)
	 	  
	  // 4) Insert into logN at the right places
	  k1 = 0;
	  for(int i = 0; i < notCondOn.size(); ++i){	    
	    if(notCondOn(i) == 1){
	      int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
	      // Age index = age - minAge
	      int a = i % nages;
	      int ageOffset = sam.confSets(s).minAge - minAgeAll;
	      int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	      Type logRec = logN.col(s)(0,y);
	      logN.col(s)(a - ageOffset,y) = simRes(k1++);
	      if(a - ageOffset == 0 && sam.confSets(s).simKeepRec){// && isForecast)
		logN.col(s)(a - ageOffset,y) = logRec;
	      }
	    }
	  }
	  // 5) Update if recruitment age is 0
	  // Need to rerun all of this:
	  bool hasZeroRec = false;
	  k1 = 0; k2 = 0;
	  for(int i = 0; i < notCondOn.size(); ++i){
	    int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
	    // Age index = age - minAge
	    int a = i % nages;
	    int ageOffset = sam.confSets(s).minAge - minAgeAll;
	    int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	    if(sam.confSets(s).minAge == 0 &&
	       sam.forecastSets(s).recModel(CppAD::Integer(sam.forecastSets(s).forecastYear(y))-1) == sam.forecastSets(s).asRecModel ){
	      hasZeroRec = true;
	      array<Type> logNa = getArray(logN, s);
	      array<Type> logFa = getArray(logF, s);
	      dataSet<Type> ds = sam.dataSets(s);
	      confSet cs = sam.confSets(s);
	      paraSet<Type> ps = paraSets(s);
	      // Overwrite (only) recruitment
	      if(a - ageOffset == 0)
		predN(i) = predNFun(ds, cs, ps, logNa, logFa, recruits(s), mortalities(s), (int)y)(0);
	      if(notCondOn(i) == 0){
		if(keepN(i) == 0){
		  avec(k1++) = 0.0;
		}else{		  
		  avec(k1++) = logN.col(s)(a - ageOffset,y) - predN(i);
		}
	      }else{
		muNew(k2++) = predN(i);
	      }
	    }
	  }
	  // Update simres and insert
	  if(hasZeroRec){
	    simRes = muNew + (vector<Type>)(meanCorrectionMat * avec) + noiseN;
	    k1 = 0;
	    for(int i = 0; i < notCondOn.size(); ++i){
	      if(notCondOn(i) == 1){
		int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
		// Age index = age - minAge
		int a = i % nages;
		int ageOffset = sam.confSets(s).minAge - minAgeAll;
		int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
		logN.col(s)(a - ageOffset,y) = simRes(k1++);
		Type logRec = logN.col(s)(0,y);
		if(a - ageOffset == 0 && sam.confSets(s).simKeepRec && isForecast){
		  logN.col(s)(a - ageOffset,y) = logRec;
		}
	      }
	    }
	  }	    
	} // End if(nCond==0){ }else if(nNotCond > 0){  }	  
      }// End doSim
      PutRNGstate();
    } // End simulate
   } // End year loop


  ///////////////////////////////////
  ////////// OBSERVATIONS //////////
  /////////////////////////////////


   // REPORT F at this point as well to get simulated forecast for MSE
   SIMULATE{
     GetRNGstate();
     for(int s = 0; s < logF.cols(); ++s){
       oftmp<Type> of(this->do_simulate);
       matrix<Type> logFs = logF.col(s);
       REPORT_F(logFs,(&of));
       ADREPORT_F(logFs,(&of));
       ofAll.addToReport(of.report,s);
       moveADREPORT(&of,this,s);
     }
     for(int s = 0; s < nStocks; ++s){
       oftmp<Type> of(this->do_simulate);     
       array<Type> logFa = getArray(logF, s);
       array<Type> lfs = getArray(logitFseason,s);
       mortalities(s) = MortalitySet<Type>(sam.dataSets(s), sam.confSets(s), paraSets(s), logFa, lfs);
       //reportMort((MortalitySet<Type>)mortalities(s),&of);
       MortalitySet<Type> mort = mortalities(s);
       REPORT_F(mort, (&of));
       ofAll.addToReport(of.report,s);
       moveADREPORT(&of,this,s);
     }
     GetRNGstate();
   }


  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    array<Type> logPa = getArray(logP, s);
    array<Type> logitFSa = getArray(logitFseason, s);
    data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);

    Type tmp = nllObs(ds, cs, ps, sam.forecastSets(s), logNa, logFa, logPa,logitFSa, recruits(s), mortalities(s), keepTmp, reportingLevel,  &of);
    if(!skip_stock_observations)
      ans += tmp;
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }

  
  ans += sharedObservation(sharedObs,
  			   sam.dataSets,
  			   sam.confSets,
  			   paraSets,
  			   sam.forecastSets,
  			   logF,
  			   logN,
  			   logP,
			   logitFseason,
  			   mortalities,
  			   // shared_logSdObs,
			   stockAreas,
			   Parea,
			   logContamination,
			   keyContamination,
  			   minYearAll,
  			   minAgeAll,
  			   keep(keep.size()-1),
  			   this);


  ///////////////////////////////
  ////////// GENETICS //////////
  /////////////////////////////

  ans += nllGenetics(sharedObs,
  		     sam.dataSets,
  		     sam.confSets,
  		     paraSets,
  		     sam.forecastSets,
  		     logF,
  		     logN,
  		     logP,
		     logitFseason,
  		     mortalities,
  		     genParSet,
  		     geneticsData,
  		     logGst,
  		     logGtrip,
		     stockAreas,
		     Parea,
		     logContamination,
		     keyContamination,
  		     maxAgeAll,
  		     minAgeAll,
		     this);

  ans += nllContamination(logContamination,
			  keyContamination,
			  meanContamination,
			  transfRhoContamination,
			  logSdContamination);


  ///////////////////////////////////////
  ////////// REFERENCE POINTS //////////
  /////////////////////////////////////

  for(int s = 0; s < nStocks; ++s){
    oftmp<Type> of(this->do_simulate);
    array<Type> logNa = getArray(logN, s);
    array<Type> logFa = getArray(logF, s);
    data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);
    reportReferencePoints(ds, cs, ps, logNa, logFa, recruits(s), sam.referencepointLists(s), &of);
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }


  ////////////////////////////////////////////////////////
  ////// ADD REPORTED OBJECTS FROM stockassessment //////
  //////////////////////////////////////////////////////

  if(isDouble<Type>::value && this->current_parallel_region<0)
    moveREPORT(this->report,ofAll.report);

   ans += getTotals(sam.dataSets,
	    sam.confSets,
	    paraSets,
	    logF,
	    logN,
	    mortalities,
	    sharedObs,
	    minYearAll,
	    maxYearAll,
	    minAgeAll,
	    maxAgeAll,
	    mohn,
	    this);

  SIMULATE{
    PutRNGstate();
  }


  return ans;
}
