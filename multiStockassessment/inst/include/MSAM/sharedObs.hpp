SAM_DEPENDS(predobs)
SAM_DEPENDS(obs)
SAM_DEPENDS(mvmix)
SAM_DEPENDS(forecast)
SAM_DEPENDS(derived)
MSAM_DEPENDS(define)
MSAM_DEPENDS(convenience)
MSAM_DEPENDS(param_types)

template<class Type>
vector<Type> predOneObsPerStock(int fleet,	// obs.aux(i,1)
				int fleetType,	// obs.fleetTypes(f-1)
				int age,	// obs.aux(i,2)-confA(s).minAge
				int year, // obs.aux(i,0)
				vector<int> noYearsLAI,
				vector<dataSet<Type> >& datA,
				vector<confSet>& confA,
				vector<paraSet<Type>>& parA,
				vector<forecastSet<Type>>& forecastA,
				cmoe_matrix<Type>& logF,
				cmoe_matrix<Type>& logN,
				vector<array<Type> >& logPs,
				vector<vector<Type> >& varAlphaSCB,
				vector<MortalitySet<Type> >& mortalities,
				vector<Type> keyFleetStock)SOURCE({
  if(fleetType > 2)
    Rf_error("fleetType not implemented for shared observations");
  vector<Type> pred(datA.size());
  pred.setConstant(R_NegInf);
  for(int s = 0; s < datA.size(); ++s){
    // Skip if age or year is not used for this stock
    int y = year -  CppAD::Integer(datA(s).years(0));
    int a = age - confA(s).minAge;
    if(y < 0 || y >= datA(s).noYears + forecastA(s).nYears)
      continue;
    if(a < 0 || a > (confA(s).maxAge - confA(s).minAge))
      continue;

    array<Type> logNa = getArray(logN,s);
    array<Type> logFa = getArray(logF,s);
    array<Type> logPa = logPs(s);

    if(fleetType == 3)
      Rf_error("Fleet type 3 not implemented for shared observations yet.");
    pred(s) = predOneObs(fleet,	// obs.aux(i,1)
			 fleetType,	// obs.fleetTypes(f-1)
			 age,	// obs.aux(i,2)-confA(s).minAge
			 year, // obs.aux(i,0)
			 CppAD::Integer(datA(s).years(0)), // minYear
			 noYearsLAI(s),
			 datA(s),
			 confA(s),
			 parA(s),
			 logFa,
			 logNa,
			 logPa,
			 varAlphaSCB(s),
			 mortalities(s),
			 (Type)0.0,
			 (Type)0.0,
			 (Type)0.0,
			 (Type)0.0,
			 (Type)0.0,
			 (Type)0.0,	     // dat.aux(i,5)
			 (Type)0.0,	     // dat.aux(i,6)
			 (Type)0.0 // releaseSurvivalVec(i)
			 );
   if(keyFleetStock(s) == 0){
     pred(s) = R_NegInf;
   }else{
     pred(s) += log(keyFleetStock(s));
   }
  }
  return pred;   
				  })


MSAM_SPECIALIZATION(vector<double> predOneObsPerStock(int, int, int, int, vector<int>, vector<dataSet<double> >&, vector<confSet>&, vector<paraSet<double> >&, vector<forecastSet<double> >&, cmoe_matrix<double>&, cmoe_matrix<double>&, vector<array<double> >&, vector<vector<double> >&, vector<MortalitySet<double> >&, vector<double>));
MSAM_SPECIALIZATION(vector<TMBad::ad_aug> predOneObsPerStock(int, int, int, int, vector<int>, vector<dataSet<TMBad::ad_aug> >&, vector<confSet>&, vector<paraSet<TMBad::ad_aug> >&, vector<forecastSet<TMBad::ad_aug> >&, cmoe_matrix<TMBad::ad_aug>&, cmoe_matrix<TMBad::ad_aug>&, vector<array<TMBad::ad_aug> >&, vector<vector<TMBad::ad_aug> >&, vector<MortalitySet<TMBad::ad_aug> >&, vector<TMBad::ad_aug>));



template<class Type>
matrix<Type> predObsPerStock(shared_obs<Type>& obs,
			     vector<dataSet<Type> >& datA,
			     vector<confSet>& confA,
			     vector<paraSet<Type> >& parA,
			     vector<forecastSet<Type> >& forecastA,
			     cmoe_matrix<Type>& logF,
			     cmoe_matrix<Type>& logN,
			     vector<array<Type> >& logP,
			     vector<vector<Type> >& varAlphaSCB,
			     vector<MortalitySet<Type> >& mortalities,
			     // vector<Type> logSdObs,
			     int minYearAll,
			     int minAgeAll,
			     vector<int> noYearsLAI,
			     objective_function<Type> *of)SOURCE({
  
  matrix<Type> predPerStock(obs.aux.dim[0], datA.size());
  predPerStock.setConstant(R_NegInf);

  //   Type logzz=Type(R_NegInf);
  for(int i=0;i<obs.aux.dim[0];i++){
    
    predPerStock.row(i) = predOneObsPerStock(obs.aux(i,1), // fleet
					     obs.fleetTypes(obs.aux(i,1)-1), // fleetType
					     obs.aux(i,2), // age
					     obs.aux(i,0), // year
					     noYearsLAI,
					     datA,
					     confA,
					     parA,
					     forecastA,
					     logF,
					     logN,
					     logP,
					     varAlphaSCB,
					     mortalities,
					     (vector<Type>)obs.keyFleetStock.row(obs.aux(i,1)-1));
  }

  return predPerStock;
}
			       );



MSAM_SPECIALIZATION(matrix<double> predObsPerStock(shared_obs<double>&, vector<dataSet<double> >&, vector<confSet>&, vector<paraSet<double> >&, vector<forecastSet<double> >&, cmoe_matrix<double>&, cmoe_matrix<double>&, vector<array<double> >&, vector<vector<double> >&, vector<MortalitySet<double> >&, int, int, vector<int>, objective_function<double>*));
MSAM_SPECIALIZATION(matrix<TMBad::ad_aug> predObsPerStock(shared_obs<TMBad::ad_aug>&, vector<dataSet<TMBad::ad_aug> >&, vector<confSet>&, vector<paraSet<TMBad::ad_aug> >&, vector<forecastSet<TMBad::ad_aug> >&, cmoe_matrix<TMBad::ad_aug>&, cmoe_matrix<TMBad::ad_aug>&, vector<array<TMBad::ad_aug> >&, vector<vector<TMBad::ad_aug> >&, vector<MortalitySet<TMBad::ad_aug> >&, int, int, vector<int>, objective_function<TMBad::ad_aug>*));




template<class Type>
Type sharedObservation(shared_obs<Type>& obs,
		       vector<dataSet<Type> >& datA,
		       vector<confSet>& confA,
		       vector<paraSet<Type> >& parA,
		       vector<forecastSet<Type> >& forecastA,
		       cmoe_matrix<Type>& logF,
		       cmoe_matrix<Type>& logN,
		       cmoe_matrix<Type>& logP,
		       vector<MortalitySet<Type> >& mortalities,
		       // vector<Type> logSdObs,
		       int minYearAll,
		       int minAgeAll,
		       data_indicator<vector<Type>,Type>& keep,
		       objective_function<Type> *of)SOURCE({
  if(!obs.hasSharedObs){
    return Type(0.0);
  }


  vector<array<Type> > comps(datA.size());
  vector<vector<Type> > weekContrib(datA.size());
  vector<int> noYearsLAI(datA.size());
  for(int s = 0; s < datA.size(); ++s){
    array<Type> logPa = getArray(logP,s);
    comps(s) = scalePFun(confA(s), datA(s), logPa);
    weekContrib(s) = scaleWeekFun(parA(s), datA(s), logPa);
    noYearsLAI(s) = yearsPFun(confA(s),datA(s));
  }
      
  matrix<Type> predPerStock = predObsPerStock(obs, datA, confA, parA, forecastA, logF, logN, comps, weekContrib, mortalities, minYearAll, minAgeAll, noYearsLAI, of);
  vector<vector< MVMIX_t<Type> > > nllVecPerStock(datA.size());
  for(int s = 0; s < nllVecPerStock.size(); ++s)
    nllVecPerStock(s) = getnllVec(datA(s), confA(s), parA(s), of);

  Type nll = 0.0;

  for(int y=0;y<obs.noYears;y++){
    int totalParKey = 0;
    for(int f=0;f<obs.noFleets;f++){
      if(!((obs.fleetTypes(f)==5)||(obs.fleetTypes(f)==3))){
        if(!isNAINT(obs.idx1(f,y))){ 
	  int idxfrom=obs.idx1(f,y);
          int idxlength=obs.idx2(f,y)-idxfrom+1;
	  // Get observation vector
	  matrix<Type> predUse = predPerStock.block(idxfrom,0,idxlength,predPerStock.cols());
	  
	  // Combine to total
	  vector<Type> predObsSegment(predUse.rows());
	  predObsSegment.setConstant(R_NegInf);
	  for(int a = 0; a < predUse.rows(); ++a){ // loop over ages
	    for(int s = 0; s < predUse.cols(); ++s){ // loop over stocks
	      predObsSegment(a) = logspace_add_SAM(predObsSegment(a), predUse(a,s));
	    }
	  }

	  MVMIX_t<Type> thisNll;
	  // obs.covCombineType
	  if(obs.covCombine(f) == shared_obs<Type>::CC_Delta_LogSumExp ||
	     obs.covCombine(f) == shared_obs<Type>::CC_Delta_SumExp){
	    matrix<Type> logggAgeStock = predUse.transpose(); // (obsUse.cols(), obsUse.rows());
	    for(int a = 0; a < logggAgeStock.cols(); ++a){
	      Type ls = logspace_sum((vector<Type>)logggAgeStock.col(a));
	      vector<Type> tmp = (vector<Type>)logggAgeStock.col(a);
	      if(obs.covCombine(f) == shared_obs<Type>::CC_Delta_LogSumExp)
		tmp -= ls;
	      logggAgeStock.col(a) = tmp;
	    }
	    // // 2) Combine covariances
	    matrix<Type> loggg(logggAgeStock.cols(), logggAgeStock.size()); // Ages x (ages*stocks)
	    loggg.setZero();
	    for(int a = 0; a < logggAgeStock.cols(); ++a)
	      for(int s = 0; s < logggAgeStock.rows(); ++s){
		loggg(a, s*logggAgeStock.cols() + a) = exp(logggAgeStock(s,a));
	      }
	    matrix<Type> SigmaAll(idxlength * predPerStock.cols(),
				  idxlength * predPerStock.cols()); // (ages*stocks) x (ages*stocks)
	    SigmaAll.setZero();
	    for(int s = 0; s < nllVecPerStock.size(); ++s){
	      matrix<Type> tmpMat = nllVecPerStock(s)(f).cov();
	      SigmaAll.block(s * loggg.rows(),s * loggg.rows(), loggg.rows(), loggg.rows()) = tmpMat;
	    }
	    matrix<Type> SigmaCombine = loggg * SigmaAll * loggg.transpose();
	    // 3) Combine t probabilities
	    Type pCombine = 0;
	    for(int s = 0; s < nllVecPerStock.size(); ++s)
	      pCombine +=  confA(s).fracMixObs(f) / (Type)nllVecPerStock.size();
	    // 4) Create MVMIX
	    thisNll = MVMIX_t<Type>(SigmaCombine, pCombine);
	  }else if(obs.covCombine(f) == shared_obs<Type>::CC_Average){
	    matrix<Type> SigmaCombine = nllVecPerStock(0)(f).cov() * 0;
	    Type pCombine = 0;
	    Type d = 0;
	    for(int s = 0; s < nllVecPerStock.size(); ++s){
	      matrix<Type> tmpMat = obs.keyFleetStock(s) * nllVecPerStock(s)(f).cov();
	      SigmaCombine += tmpMat;
	      pCombine += obs.keyFleetStock(s) * confA(s).fracMixObs(f);
	      d += obs.keyFleetStock(s);
	    }
	    SigmaCombine /= d;
	    pCombine /= d;
	    // 4) Create MVMIX
	    thisNll = MVMIX_t<Type>(SigmaCombine, pCombine);
	  }else{
	    Rf_error("Covaraince combination type not implemented");
	  }
	  
          vector<Type> currentVar=thisNll.cov().diagonal();
          vector<Type> sqrtW(currentVar.size());
	  sqrtW.setConstant(1.0);
	  switch(confA(0).obsLikelihoodFlag(f)){
	  case 0: // (LN) log-Normal distribution
              
	    nll += thisNll((obs.logobs.segment(idxfrom,idxlength)-predObsSegment), keep.segment(idxfrom,idxlength)); ///sqrtW,keep.segment(idxfrom,idxlength));
	    //     nll += (log(sqrtW)*keep.segment(idxfrom,idxlength)).sum();
	      SIMULATE_F(of){
		obs.logobs.segment(idxfrom,idxlength) = predObsSegment + (thisNll.simulate()*sqrtW);
	      }
	    break;
	  case 1: // (ALN) Additive logistic-normal proportions + log-normal total numbers
	    nll +=  thisNll(obs_fun::addLogratio((vector<Type>)obs.logobs.segment(idxfrom,idxlength))-obs_fun::addLogratio((vector<Type>)predObsSegment));
	    nll += log(obs_fun::log2proportion((vector<Type>)obs.logobs.segment(idxfrom,idxlength))).sum();
	    nll -= dnorm(log(obs_fun::log2expsum((vector<Type>)obs.logobs.segment(idxfrom,idxlength))),
	     	         log(obs_fun::log2expsum((vector<Type>)predObsSegment)),
	   	         exp(parA(0).logSdLogTotalObs(totalParKey++)),true); // SHOULD BE CHANGED TO USE DELTA METHOD
	    nll += log(obs_fun::log2expsum((vector<Type>)obs.logobs.segment(idxfrom,idxlength)));
	    nll -= log(fabs(obs_fun::jacobianDet((vector<Type>)obs.logobs.segment(idxfrom,idxlength).exp())));
            nll -= obs.logobs.segment(idxfrom,idxlength).sum();
	    SIMULATE_F(of){
	      vector<Type> logProb(idxlength);
	      logProb.setZero();
	      logProb.segment(0,idxlength-1) = obs_fun::addLogratio(((vector<Type>)predObsSegment)) + thisNll.simulate();
	      Type logDenom = obs_fun::logExpSum(logProb);
	      logProb -= logDenom;
	      Type logTotal = rnorm(log(obs_fun::log2expsum((vector<Type>)predObsSegment)),
	  			    exp(parA(0).logSdLogTotalObs(totalParKey++)));
	      obs.logobs.segment(idxfrom,idxlength) = logProb + logTotal; 
	    }
	    break;
	  default:
	    Rf_error("Unknown obsLikelihoodFlag");
	  }
        }
 
      }   
    }  
  }
  return nll; //nll;
}
			 );


MSAM_SPECIALIZATION(double sharedObservation(shared_obs<double>&, vector<dataSet<double> >&, vector<confSet>&, vector<paraSet<double> >&, vector<forecastSet<double> >&, cmoe_matrix<double>&, cmoe_matrix<double>&, cmoe_matrix<double>&, vector<MortalitySet<double> >&, int, int, data_indicator<vector<double>, double>&, objective_function<double>*));
MSAM_SPECIALIZATION(TMBad::ad_aug sharedObservation(shared_obs<TMBad::ad_aug>&, vector<dataSet<TMBad::ad_aug> >&, vector<confSet>&, vector<paraSet<TMBad::ad_aug> >&, vector<forecastSet<TMBad::ad_aug> >&, cmoe_matrix<TMBad::ad_aug>&, cmoe_matrix<TMBad::ad_aug>&, cmoe_matrix<TMBad::ad_aug>&, vector<MortalitySet<TMBad::ad_aug> >&, int, int, data_indicator<vector<TMBad::ad_aug>, TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));
