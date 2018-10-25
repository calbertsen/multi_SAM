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

  
  PARAMETER_CMOE_VECTOR(logFpar);
  PARAMETER_CMOE_VECTOR(logQpow);
  PARAMETER_CMOE_VECTOR(logSdLogFsta);
  PARAMETER_CMOE_VECTOR(logSdLogN);
  PARAMETER_CMOE_VECTOR(logSdLogObs);
  PARAMETER_CMOE_VECTOR(logSdLogTotalObs);
  PARAMETER_CMOE_VECTOR(transfIRARdist);
  PARAMETER_CMOE_VECTOR(sigmaObsParUS);
  PARAMETER_CMOE_VECTOR(rec_loga);
  PARAMETER_CMOE_VECTOR(rec_logb);
  PARAMETER_CMOE_VECTOR(itrans_rho);
  PARAMETER_CMOE_VECTOR(logScale);
  PARAMETER_CMOE_VECTOR(logitReleaseSurvival);
  PARAMETER_CMOE_VECTOR(logitRecapturePhi);

  PARAMETER_CMOE_MATRIX(logF); 
  PARAMETER_CMOE_MATRIX(logN);
  PARAMETER_CMOE_VECTOR(missing);

  ////////////////////////////////////////
  ////////// Multi SAM specific //////////
  ////////////////////////////////////////
  PARAMETER_VECTOR(RE)

    
    int nStocks = logF.cols();
  vector<paraSet<Type> > paraSets(nStocks);

  for(int s = 0; s < nStocks; ++s){
    paraSets(s).logFpar = logFpar.col(s);  
    paraSets(s).logQpow = logQpow.col(s);  
    paraSets(s).logSdLogFsta = logSdLogFsta.col(s);
    paraSets(s).logSdLogN = logSdLogN.col(s);
    paraSets(s).logSdLogObs = logSdLogObs.col(s);
    paraSets(s).logSdLogTotalObs = logSdLogTotalObs.col(s);
    paraSets(s).transfIRARdist = transfIRARdist.col(s); //transformed distances for IRAR cor obs structure
    paraSets(s).sigmaObsParUS = sigmaObsParUS.col(s); //choleski elements for unstructured cor obs structure
    paraSets(s).rec_loga = rec_loga.col(s);
    paraSets(s).rec_logb = rec_logb.col(s);
    paraSets(s).itrans_rho = itrans_rho.col(s);  
    paraSets(s).logScale = logScale.col(s);
    paraSets(s).logitReleaseSurvival = logitReleaseSurvival.col(s);
    paraSets(s).logitRecapturePhi = logitRecapturePhi.col(s);    
  }

  // patch missing 
  int idxmis=0;
  for(int s = 0; s < nStocks; ++s)
    for(int i = 0; i < sam.dataSets(s).nobs; ++i){
      if(isNA(sam.dataSets(s).logobs(i))){
  	sam.dataSets(s).logobs(i)=missing(idxmis++);
      }    
    }
  
  // vector<vector<Type> > R(nStocks);
  // vector<vector<Type> > logR(nStocks);
  // vector<vector<Type> > ssb(nStocks);
  // vector<vector<Type> > logssb(nStocks);
  // vector<vector<Type> > fbar(nStocks);
  // vector<vector<Type> > logfbar(nStocks);
  // vector<vector<Type> > cat(nStocks);
  // vector<vector<Type> > logCatch(nStocks);
  // vector<vector<Type> > varLogCatch(nStocks);
  // vector<vector<Type> > fsb(nStocks);
  // vector<vector<Type> > logfsb(nStocks);
  // vector<vector<Type> > tsb(nStocks);
  // vector<vector<Type> > logtsb(nStocks);
  // vector<vector<Type> > predObs(nStocks);


  // // Reported inside functions
  // vector<vector<matrix<Type> > > obsCov(nStocks);
  
  // for(int s = 0; s < nStocks; ++s){
  //   array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
  //   logNa = logN.col(s);
  //   array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
  //   logFa = logF.col(s);
  //   R(s) = rFun(logNa);
  //   logR(s) = log(R(s));  
  //   ssb(s) = ssbFun(sam.dataSets(s), sam.confSets(s), logNa, logFa);
  //   logssb(s) = log(ssb(s));
  //   fbar(s) = fbarFun(sam.confSets(s), logFa);
  //   logfbar(s) = log(fbar(s));
  //   cat(s) = catchFun(sam.dataSets(s), sam.confSets(s), logNa, logFa);
  //   logCatch(s) = log(cat(s));
  //   varLogCatch(s) = varLogCatchFun(sam.dataSets(s), sam.confSets(s), logNa, logFa, paraSets(s));
  //   fsb(s) = fsbFun(sam.dataSets(s), sam.confSets(s), logNa, logFa);
  //   logfsb(s) = log(fsb(s));
  //   tsb(s) = tsbFun(sam.dataSets(s), sam.confSets(s), logNa);
  //   logtsb(s) = log(tsb(s));
  //   predObs(s) = predObsFun(sam.dataSets(s), sam.confSets(s), paraSets(s), logNa, logFa, logssb(s), logfsb(s), logCatch(s));
  // }


  

  Type ans = 0; //negative log-likelihood
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

    data_indicator<vector<Type>,Type> keep(sam.dataSets(s).logobs);
    ans += nllF(sam.confSets(s), paraSets(s), logFa, keep, &of);
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
  matrix<Type> ncov;

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
      A(i,i) = Type(0.0);
    }
  }
  REPORT(A);

  matrix<Type> L = constructL(nages,nAreas,RE,cons);
    
  // Get covariance matrix with arbitrary scale
  matrix<Type> Sigma = L * L.transpose();
  Sigma = cov2cor(Sigma);
  matrix<Type> SigmaTmp = Sigma;
    
  if(usePartialCor){
    Sigma = cov2cor((matrix<Type>)Sigma.inverse()); //pcor2cor(Sigma);
  }

  ncov = (matrix<Type>)(A * Sigma) * A;
  REPORT(Sigma);
  REPORT(SigmaTmp);
  REPORT(L);
  REPORT(ncov)

    ////////////////////////////////////
    ////// CALCULATE CONTRIBUTION //////
    ////////////////////////////////////

    density::MVNORM_t<Type> neg_log_densityN(ncov);
  // Loop over time
  for(int yall = 0; yall < maxYearAll - minYearAll + 1; ++yall){

    // Vector for predictions
    vector<Type> predN(ncov.rows());
    predN.setZero();
    vector<Type> newN(ncov.rows());
    newN.setZero();
    // Determine which ages to use with keep vector
    vector<Type> keep(ncov.rows());
    keep.setZero();
    // Loop over stocks
    for(int s = 0; s < nAreas; ++s){
      int ageOffset = sam.confSets(s).minAge - minAgeAll;
      int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
      
      array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
      logNa = logN.col(s);
      array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
      logFa = logF.col(s);
      if(y > 0 && y < sam.dataSets(s).noYears){
	vector<Type> predNnz = predNFun(sam.dataSets(s), sam.confSets(s), paraSets(s), logNa, logFa, y);
	keep.segment(s * nages + ageOffset,predNnz.size()) = 1.0;
	predN.segment(s * nages + ageOffset,predNnz.size()) = predNnz;
	newN.segment(s * nages + ageOffset,predNnz.size()) = logNa.col(y);
      }
      REPORT(keep);REPORT(predN);REPORT(newN);
    }
    if(keep.sum()>0)
      ans+=neg_log_densityN(newN-predN, keep);
    SIMULATE{
      // R_NaReal
      // 1) Simulate all values
      vector<Type> Ntmp = predN + neg_log_densityN.simulate();
      // 2) Insert in logN if conf.simFlag==0
      for(int i = 0; i < keep.size(); ++i){
	if(keep(i) > 0){
	  int s = (int)i / (int)nages; // must be integer division: A.rows() = nages * nAreas
	  // Age index = age - minAge
	  int a = i % nages;
	  int y = yall - CppAD::Integer(sam.dataSets(s).years(0) - minYearAll);
	  if(sam.confSets(s).simFlag == 0)
	    logN.col(s)(a,y) = Ntmp(i);
	}
      } // End keep loop
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
    data_indicator<vector<Type>,Type> keep(sam.dataSets(s).logobs);
    ans += nllObs(sam.dataSets(s), sam.confSets(s), paraSets(s), logNa, logFa, keep,  &of);
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
