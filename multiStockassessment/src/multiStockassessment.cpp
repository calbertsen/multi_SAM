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

#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include <TMB.hpp>
#include <SAM.hpp>
#include "../inst/include/multiSAM.hpp"


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
  using CppAD::abs;


  DATA_STRUCT(sam,sam_data);

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


  // Related to residuals
  DATA_VECTOR(fake_obs);
  DATA_IVECTOR(fake_stock);
  DATA_IVECTOR(fake_indx);
  DATA_VECTOR_INDICATOR(fake_keep, fake_obs);
  DATA_INTEGER(doResiduals);

  // Create a keep indicator for each stock
  vector<data_indicator<vector<Type>,Type> > keep(sam.dataSets.size());
  
  for(int s = 0; s < sam.dataSets.size(); ++s){
    data_indicator<vector<Type>,Type> keepTmp(sam.dataSets(s).logobs, true);
    keep(s) = keepTmp;
  }
  
  if(doResiduals){
    for(int i = 0; i < fake_obs.size(); ++i){
      // Overwrite true observations with residual observations
      sam.dataSets(fake_stock(i)).logobs(fake_indx(i)) = fake_obs(i);
      // Overwrite keep value
      keep(fake_stock(i))(fake_indx(i)) = fake_keep(i);
    }
  }


  
  
  PARAMETER_CMOE_VECTOR(logFpar);
  PARAMETER_CMOE_VECTOR(logQpow);
  PARAMETER_CMOE_VECTOR(logSdLogFsta);
  PARAMETER_CMOE_VECTOR(logSdLogN);
  PARAMETER_CMOE_VECTOR(logSdLogObs);
  PARAMETER_CMOE_VECTOR(logSdLogTotalObs);
  PARAMETER_CMOE_VECTOR(transfIRARdist);
  PARAMETER_CMOE_VECTOR(sigmaObsParUS);
  PARAMETER_CMOE_VECTOR(rec_pars);
  PARAMETER_CMOE_VECTOR(itrans_rho);
  PARAMETER_CMOE_VECTOR(logScale);
  PARAMETER_CMOE_VECTOR(logitReleaseSurvival);
  PARAMETER_CMOE_VECTOR(logitRecapturePhi);
  PARAMETER_CMOE_VECTOR(sepFalpha);
  PARAMETER_CMOE_VECTOR(sepFlogitRho);
  PARAMETER_CMOE_VECTOR(sepFlogSd);
  PARAMETER_CMOE_VECTOR(predVarObs);

  PARAMETER_CMOE_VECTOR(sepFalpha);
  PARAMETER_CMOE_VECTOR(sepFlogitRho);
  PARAMETER_CMOE_VECTOR(sepFlogSd);

  // Forecast FMSY
  PARAMETER_VECTOR(logFScaleMSY);
  PARAMETER_VECTOR(implicitFunctionDelta);

  // YPR reference points
  PARAMETER_VECTOR(logScaleFmsy);
  PARAMETER_VECTOR(logScaleFmax);
  PARAMETER_VECTOR(logScaleF01);
  PARAMETER_VECTOR(logScaleFcrash);
  PARAMETER_VECTOR(logScaleFext);
  PARAMETER_CMOE_VECTOR(logScaleFxPercent);
  PARAMETER_VECTOR(logScaleFlim);


  PARAMETER_CMOE_MATRIX(logF); 
  PARAMETER_CMOE_MATRIX(logN);
  PARAMETER_CMOE_VECTOR(missing);

  
  
  ////////////////////////////////////////
  ////////// Multi SAM specific //////////
  ////////////////////////////////////////
  PARAMETER_VECTOR(RE);

  int nStocks = logF.cols();
  
  vector<paraSet<Type> > paraSets(nStocks);

  for(int s = 0; s < nStocks; ++s){
    paraSets(s).logFpar = logFpar.col(s);  
    paraSets(s).logQpow = logQpow.col(s);  
    paraSets(s).logSdLogFsta = logSdLogFsta.col(s);
    paraSets(s).logSdLogN = logSdLogN.col(s);
    paraSets(s).logSdLogObs = logSdLogObs.col(s);
    paraSets(s).logSdLogTotalObs = logSdLogTotalObs.col(s);
    paraSets(s).transfIRARdist = transfIRARdist.col(s);
    paraSets(s).sigmaObsParUS = sigmaObsParUS.col(s);
    paraSets(s).rec_pars = rec_pars.col(s);
    paraSets(s).itrans_rho = itrans_rho.col(s);  
    paraSets(s).logScale = logScale.col(s);
    paraSets(s).logitReleaseSurvival = logitReleaseSurvival.col(s);
    paraSets(s).logitRecapturePhi = logitRecapturePhi.col(s);
    paraSets(s).sepFalpha = sepFalpha.col(s);
    paraSets(s).sepFlogitRho = sepFlogitRho.col(s);
    paraSets(s).sepFlogSd = sepFlogSd.col(s);
    paraSets(s).predVarObs = predVarObs.col(s);    

    // Forecast FMSY
    paraSets(s).logFScaleMSY = logFScaleMSY(s);
    paraSets(s).implicitFunctionDelta = implicitFunctionDelta(s);
    // YPR reference points
    paraSets(s).logScaleFmsy = logScaleFmsy(s);
    paraSets(s).logScaleFmax = logScaleFmax(s);
    paraSets(s).logScaleF01 = logScaleF01(s);
    paraSets(s).logScaleFcrash = logScaleFcrash(s);
    paraSets(s).logScaleFext = logScaleFext(s);
    paraSets(s).logScaleFxPercent = logScaleFxPercent.col(s);
    paraSets(s).logScaleFlim = logScaleFlim(s);
  }

  
  ////////////////////////////////////////
  /////////// Prepare forecast ///////////
  ////////////////////////////////////////
  
  for(int s = 0; s < nStocks; ++s){
    // Resize arrays
    prepareForForecast(sam.dataSets(s));
    // Calculate forecast
    array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
    logNa = logN.col(s);
    array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
    logFa = logF.col(s);
    sam.dataSets(s).forecast.calculateForecast(logFa,logNa, sam.dataSets(s), sam.confSets(s), paraSets(s));
  }


  Type ans = 0; //negative log-likelihood
  
  // patch missing 
  int idxmis = 0;
  for(int s = 0; s < nStocks; ++s)
    for(int i = 0; i < sam.dataSets(s).nobs; ++i){
      if(isNA(sam.dataSets(s).logobs(i))){
  	sam.dataSets(s).logobs(i) = missing(idxmis++);
      }    
    }
  
  // add wide prior for missing random effects, but _only_ when computing ooa residuals
  if(doResiduals){
    Type huge = 10.0;
    for (int i = 0; i < missing.size(); i++)
      ans -= dnorm(missing(i), Type(0.0), huge, true);  
  } 
  
  ofall<Type> ofAll(nStocks);
  ////////////////////////////////
  ////////// F PROCESS //////////
  //////////////////////////////
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
    ans += nllF(ds, cs, ps, logFa, keepTmp, &of);
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
    // If simulate -> move grab new logF values and move them to the right place!         
    logF.col(s) = logFa.matrix();
  }

  ////////////////////////////////
  ////////// N PROCESS //////////
  //////////////////////////////
  int nAreas = sam.dataSets.size();

  int nages = (maxAgeAll - minAgeAll + 1);
 
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

    // Vector for predictions
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
      
      array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
      logNa = logN.col(s);
      array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
      logFa = logF.col(s);
      if(y > 0 && y < ds.noYears + ds.forecast.nYears){
	vector<Type> predNnz = predNFun(ds, cs, ps, logNa, logFa, (int)y);
	keepN.segment(s * nages + ageOffset,predNnz.size()) = 1.0;
	predN.segment(s * nages + ageOffset,predNnz.size()) = predNnz;
	newN.segment(s * nages + ageOffset,predNnz.size()) = logNa.col(y);
	
      }
    }
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
	   y < ds.noYears + ds.forecast.nYears &&
	   ds.forecast.nYears > 0 &&
	   ds.forecast.recModel(CppAD::Integer(ds.forecast.forecastYear(y))-1) != ds.forecast.asRecModel &&
	   ds.forecast.forecastYear(y) > 0){
	  Nscale(s * nages + ageOffset) = sqrt(ds.forecast.logRecruitmentVar) / sqrt(ncov(s * nages + ageOffset,s * nages + ageOffset));
	  predNTmp(s * nages + ageOffset) = ds.forecast.logRecruitmentMedian;
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
      for(int s = 0; s < nAreas; ++s)
	if(sam.confSets(s).simFlag(1) == 0){
	  doSim = true;
	  break;
	}
      if(doSim){
	// 2) Get conditional distribution of stocks with simFlag == 0 (and keepN = 1)
	vector<int> notCondOn(ncov.cols());
	notCondOn.setZero();
	for(int i = 0; i < notCondOn.size(); ++i){
	  int s = (int)i / (int)nages;
	  if(sam.confSets(s).simFlag(1) == 0 && keepN(i) == 1)
	    notCondOn(i) = 1;
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
    array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
    logNa = logN.col(s);
    array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
    logFa = logF.col(s);
    data_indicator<vector<Type>,Type> keepTmp = keep(s);
    dataSet<Type> ds = sam.dataSets(s);
    confSet cs = sam.confSets(s);
    paraSet<Type> ps = paraSets(s);
    ans += nllObs(ds, cs, ps, logNa, logFa, keepTmp,  &of);
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }

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
    ans += nllReferencepoints(ds, cs, ps, logNa, logFa, &of);
    ofAll.addToReport(of.report,s);
    moveADREPORT(&of,this,s);
  }

  ////////////////////////////////////////////////////////
  ////// ADD REPORTED OBJECTS FROM stockassessment //////
  //////////////////////////////////////////////////////
  if(isDouble<Type>::value && this->current_parallel_region<0)
    moveREPORT(this->report,ofAll.report);
  
  return ans;
}
