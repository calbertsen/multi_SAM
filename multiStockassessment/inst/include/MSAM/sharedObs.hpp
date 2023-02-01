SAM_DEPENDS(predobs)
SAM_DEPENDS(obs)
SAM_DEPENDS(mvmix)
SAM_DEPENDS(forecast)
SAM_DEPENDS(derived)
MSM_DEPENDS(define)
MSM_DEPENDS(convenience)
MSM_DEPENDS(param_types)

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
				cmoe_3darray<Type>& logitFseason,
				vector<vector<Type> >& varAlphaSCB,
				vector<MortalitySet<Type> >& mortalities,
				vector<Type> auxData,
				matrix<Type> Parea,
				// vector<Type> logssb,
				// vector<Type> logtsb,
				// vector<Type> logfsb,
				// vector<Type> logCatch,
				// vector<Type> logLand,
				// vectpr<Type> logfbar,
				// vector<Type> releaseSurvival, // releaseSurvivalVec(i)
				vector<Type> keyFleetStock)SOURCE({
  if(fleetType > 2 && fleetType < 80)
    Rf_error("fleetType not implemented for shared observations");
  vector<Type> pred(datA.size());
  pred.setConstant(R_NegInf);

  if(fleetType == 80){		// Total Catch/Landing proportions per season
    SAM_ASSERT(auxData.size() >= 6,"aux is not large enough for fleet type 80");
    // Calculate proportion in pred(0), rest remain neginf
    Type Ctotal = 0.0;
    Type Cseason = 0.0;
    int flt = CppAD::Integer(auxData(0))-1; // TODO: allow auxData(0)==0 to sum over all fleets
    for(int s = 0; s < datA.size(); ++s){
      int y = year - CppAD::Integer(datA(s).years(0));
      int aMin = datA(s).minAgePerFleet(flt);
      int aMax = datA(s).maxAgePerFleet(flt);
      for(int aa = aMin - confA(s).minAge; aa < aMax - confA(s).minAge; ++aa){
	//Type Cttmp = exp(logN.col(s)(aa,y)) * mortalities(s).CIF(flt,aa,y,datA(s).sampleTimesStart(flt),datA(s).sampleTimesEnd(flt));
	Type Cttmp = exp(logN.col(s)(aa,y) + mortalities(s).logFleetSurvival_before(aa,y,flt) + log(mortalities(s).fleetCumulativeIncidence(aa,y,flt)));
	Type Cstmp = exp(logN.col(s)(aa,y)) * mortalities(s).partialCIF(flt,aa,y, auxData(1), auxData(2));
	// 0: Catch numbers
	if(CppAD::Integer(auxData(3)) == 1){ // 1: Catch weight
	  Cttmp *= datA(s).catchMeanWeight(y,aa, flt);
	  Cstmp *= datA(s).catchMeanWeight(y,aa, flt);
	}else if(CppAD::Integer(auxData(3)) == 2){ // 2: Landing numbers
	  Cttmp *= datA(s).landFrac(y,aa, flt);
	  Cstmp *= datA(s).landFrac(y,aa, flt);
	}else if(CppAD::Integer(auxData(3)) == 3){ // 3: Landing weight
	  Cttmp *= datA(s).landFrac(y,aa, flt) * datA(s).landMeanWeight(y,aa, flt);
	  Cstmp *= datA(s).landFrac(y,aa, flt) * datA(s).landMeanWeight(y,aa, flt);
	}
	Ctotal += Cttmp;
	Cseason += Cstmp;
      }
    }
    pred(0) = log(Cseason) - log(Ctotal);
    
  }else if(fleetType == 81){ // Age-wise Catch/Landing proportions per season
    // Sum up catch/landings and convert to proportions later
    SAM_ASSERT(auxData.size() >= 6,"aux is not large enough for fleet type 81");

    // aux: (0) catch fleet, (1) this season start, (2) this season end, (3) catch Type (4) season number, (5) number of seasons
    // For combining stocks, relative catch weights/landing fractions are important
    int flt = CppAD::Integer(auxData(0))-1;
    Type Ctotal = 0.0;
    Type Cseason = 0.0;
    for(int s = 0; s < datA.size(); ++s){
      int y = year - CppAD::Integer(datA(s).years(0));
      int a=age-confA(s).minAge;
      // TODO: allow auxData(0)==0  to sum over all fleets
      //Type Cttmp = exp(logN.col(s)(a,y)) * mortalities(s).CIF(flt,a,y,datA(s).sampleTimesStart(flt),datA(s).sampleTimesEnd(flt));
      Type Cttmp = exp(logN.col(s)(a,y) + mortalities(s).logFleetSurvival_before(a,y,flt) + log(mortalities(s).fleetCumulativeIncidence(a,y,flt)));
      Type Cstmp = exp(logN.col(s)(a,y)) * mortalities(s).partialCIF(flt,a,y, auxData(1), auxData(2));
      // 0: Catch numbers
      if(CppAD::Integer(auxData(3)) == 1){ // 1: Catch weight
	Cttmp *= datA(s).catchMeanWeight(y,a, flt);
	Cstmp *= datA(s).catchMeanWeight(y,a, flt);
      }else if(CppAD::Integer(auxData(3)) == 2){ // 2: Landing numbers
	Cttmp *= datA(s).landFrac(y,a, flt);
	Cstmp *= datA(s).landFrac(y,a, flt);
      }else if(CppAD::Integer(auxData(3)) == 3){ // 3: Landing weight
	Cttmp *= datA(s).landFrac(y,a, flt) * datA(s).landMeanWeight(y,a, flt);
	Cstmp *= datA(s).landFrac(y,a, flt) * datA(s).landMeanWeight(y,a, flt);
      }
      Ctotal += Cttmp;
      Cseason += Cstmp;
    }
    pred(0) = log(Cseason) - log(Ctotal);
    
  }else if(fleetType == 90){	// Total Stock composition in catch/landing
    SAM_ASSERT(auxData.size() >= 5,"aux is not large enough for fleet type 90");
    // Sum up catch/landings and convert to proportions later    
    // NOTE: Observations come in one at a time, so result should be -Inf for all stocks but one?
    // - All stocks needs to be calculated to get proportion. Could postpone dividing by sum? But last stock will not be in the data
    // aux: (0) catch fleet, (1) this season start, (2) this season end, (3) catch Type (4) stock number
    vector<Type> Cstock(datA.size());
    Cstock.setZero();
    Type Ctotal = 0.0;
    int flt0 = CppAD::Integer(auxData(0))-1; // TODO: allow auxData(0)==0 to sum over all fleets
    int flt1 = flt0;
    if(flt0 < 0){
      flt0 = 0;
      flt1 = datA(0).fleetTypes.size()-1;
    }
    for(int s = 0; s < datA.size(); ++s){      
      int y = year - CppAD::Integer(datA(s).years(0));
      // Calculate
      for(int flt = flt0; flt <= flt1; ++flt){
	if(datA(s).fleetTypes(flt) == 0){
	  int aMin = datA(s).minAgePerFleet(flt);
	  int aMax = datA(s).maxAgePerFleet(flt);
	  for(int aa = aMin - confA(s).minAge; aa < aMax - confA(s).minAge; ++aa){
	    Type Cstmp = exp(logN.col(s)(aa,y)) * mortalities(s).partialCIF(flt,aa,y, auxData(1), auxData(2));
	    // 0: Catch numbers
	    if(CppAD::Integer(auxData(3)) == 1){ // 1: Catch weight
	      Cstmp *= datA(s).catchMeanWeight(y,aa, flt);
	    }else if(CppAD::Integer(auxData(3)) == 2){ // 2: Landing numbers
	      Cstmp *= datA(s).landFrac(y,aa, flt);
	    }else if(CppAD::Integer(auxData(3)) == 3){ // 3: Landing weight
	      Cstmp *= datA(s).landFrac(y,aa, flt) * datA(s).landMeanWeight(y,aa, flt);
	    }
	    Ctotal += Cstmp;
	    Cstock(s) += Cstmp;
	  }
	}
      }
    }
    for(int s = 0; s < datA.size(); ++s){

      pred(s) = log(squeeze(Cstock(s) / Ctotal));
    }

  }else if(fleetType == 92){
    SAM_ASSERT(auxData.size() >= 5,"aux is not large enough for fleet type 90");

    Type Ctotal = 0.0;
    Type Carea = 0.0;
    int flt = CppAD::Integer(auxData(0))-1; // TODO: allow auxData(0)==0 to sum over all fleets
    for(int s = 0; s < datA.size(); ++s){
      int y = year - CppAD::Integer(datA(s).years(0));
      int aMin = datA(s).minAgePerFleet(flt);
      int aMax = datA(s).maxAgePerFleet(flt);
      for(int aa = aMin - confA(s).minAge; aa < aMax - confA(s).minAge; ++aa){
	//Type Cttmp = exp(logN.col(s)(aa,y)) * mortalities(s).CIF(flt,aa,y,datA(s).sampleTimesStart(flt),datA(s).sampleTimesEnd(flt));
	// Type Cttmp = exp(logN.col(s)(aa,y) + mortalities(s).logFleetSurvival_before(aa,y,flt) + log(mortalities(s).fleetCumulativeIncidence(aa,y,flt)));
	  Type Cttmp = exp(logN.col(s)(aa,y)) * mortalities(s).partialCIF(flt,aa,y, auxData(1), auxData(2));
	  // 0: Catch numbers
	  if(CppAD::Integer(auxData(3)) == 1){ // 1: Catch weight
	    Cttmp *= datA(s).catchMeanWeight(y,aa, flt);
	  }else if(CppAD::Integer(auxData(3)) == 2){ // 2: Landing numbers
	    Cttmp *= datA(s).landFrac(y,aa, flt);
	  }else if(CppAD::Integer(auxData(3)) == 3){ // 3: Landing weight
	    Cttmp *= datA(s).landFrac(y,aa, flt) * datA(s).landMeanWeight(y,aa, flt);
	  }
	  Carea += Parea(CppAD::Integer(auxData(4)),s) * Cttmp;
	  Ctotal += Cttmp;
	}
    }
    pred(0) = log(Carea) - log(Ctotal);
    
    
  }else{ // if(fleetType <= 2){
    vector<Type> NAs(datA.size());
    NAs.setConstant(R_NaReal);
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
      array<Type> logitFSa = getArray(logitFseason,s);
      array<Type> logPa = logPs(s);
    
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
			   logitFSa,
			   varAlphaSCB(s),
			   mortalities(s),
			   (Type)NAs(s),//logssb(s),
			   (Type)NAs(s),//logtsb(s),
			   (Type)NAs(s),//logfsb(s),
			   (Type)NAs(s),//logCatch(s),
			   (Type)NAs(s),//logLand(s),
			   (Type)NAs(s),//logfbar(s),
			   (Type)NAs(s),//releaseSurvival(s),
			   auxData
			   );
      if(keyFleetStock(s) == 0){
	pred(s) = R_NegInf;
      }else{
	pred(s) += log(keyFleetStock(s));
      }
    }
  // }else{
  //   Rf_error("fleetType not implemented for shared observations");
  }
  return pred;
				  })


				  MSM_SPECIALIZATION(vector<double> predOneObsPerStock(int, int, int, int, vector<int>, vector<dataSet<double> >&, vector<confSet>&, vector<paraSet<double> >&, vector<forecastSet<double> >&, cmoe_matrix<double>&, cmoe_matrix<double>&, vector<array<double> >&,cmoe_3darray<double>&, vector<vector<double> >&, vector<MortalitySet<double> >&, vector<double>,matrix<double>, vector<double>));
MSM_SPECIALIZATION(vector<TMBad::ad_aug> predOneObsPerStock(int, int, int, int, vector<int>, vector<dataSet<TMBad::ad_aug> >&, vector<confSet>&, vector<paraSet<TMBad::ad_aug> >&, vector<forecastSet<TMBad::ad_aug> >&, cmoe_matrix<TMBad::ad_aug>&, cmoe_matrix<TMBad::ad_aug>&, vector<array<TMBad::ad_aug> >&,cmoe_3darray<TMBad::ad_aug>&, vector<vector<TMBad::ad_aug> >&, vector<MortalitySet<TMBad::ad_aug> >&, vector<TMBad::ad_aug>,matrix<TMBad::ad_aug>, vector<TMBad::ad_aug>));



template<class Type>
matrix<Type> predObsPerStock(shared_obs<Type>& obs,
			     vector<dataSet<Type> >& datA,
			     vector<confSet>& confA,
			     vector<paraSet<Type> >& parA,
			     vector<forecastSet<Type> >& forecastA,
			     cmoe_matrix<Type>& logF,
			     cmoe_matrix<Type>& logN,
			     vector<array<Type> >& logP,
			     cmoe_3darray<Type>& logitFseason,
			     vector<vector<Type> >& varAlphaSCB,
			     vector<MortalitySet<Type> >& mortalities,
			     matrix<Type>& Parea,
			     // vector<Type> logSdObs,
			     int minYearAll,
			     int minAgeAll,
			     vector<int> noYearsLAI,
			     objective_function<Type> *of)SOURCE({
  
  matrix<Type> predPerStock(obs.aux.dim[0], datA.size());
  predPerStock.setConstant(R_NegInf);

  //   Type logzz=Type(R_NegInf);
  for(int i=0;i<obs.aux.dim[0];i++){
    vector<Type> auxData(obs.auxData.cols());
    for(int q = 0; q < auxData.size(); ++q)
      auxData(q) = obs.auxData(i,q);
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
					     logitFseason,
					     varAlphaSCB,
					     mortalities,
					     auxData,
					     Parea,
					     (vector<Type>)obs.keyFleetStock.row(obs.aux(i,1)-1));
  }
  REPORT_F(predPerStock,of);
  return predPerStock;
}
			       );



MSM_SPECIALIZATION(matrix<double> predObsPerStock(shared_obs<double>&, vector<dataSet<double> >&, vector<confSet>&, vector<paraSet<double> >&, vector<forecastSet<double> >&, cmoe_matrix<double>&, cmoe_matrix<double>&, vector<array<double> >&, cmoe_3darray<double>&, vector<vector<double> >&, vector<MortalitySet<double> >&,matrix<double>&, int, int, vector<int>, objective_function<double>*));
MSM_SPECIALIZATION(matrix<TMBad::ad_aug> predObsPerStock(shared_obs<TMBad::ad_aug>&, vector<dataSet<TMBad::ad_aug> >&, vector<confSet>&, vector<paraSet<TMBad::ad_aug> >&, vector<forecastSet<TMBad::ad_aug> >&, cmoe_matrix<TMBad::ad_aug>&, cmoe_matrix<TMBad::ad_aug>&, vector<array<TMBad::ad_aug> >&, cmoe_3darray<TMBad::ad_aug>&, vector<vector<TMBad::ad_aug> >&, vector<MortalitySet<TMBad::ad_aug> >&,matrix<TMBad::ad_aug>&, int, int, vector<int>, objective_function<TMBad::ad_aug>*));




template<class Type>
Type sharedObservation(shared_obs<Type>& obs,
		       vector<dataSet<Type> >& datA,
		       vector<confSet>& confA,
		       vector<paraSet<Type> >& parA,
		       vector<forecastSet<Type> >& forecastA,
		       cmoe_matrix<Type>& logF,
		       cmoe_matrix<Type>& logN,
		       cmoe_matrix<Type>& logP,
		       cmoe_3darray<Type>& logitFseason,
		       vector<MortalitySet<Type> >& mortalities,
		       matrix<Type>& Parea,
		       // vector<Type> logSdObs,
		       int minYearAll,
		       int minAgeAll,
		       data_indicator<vector<Type>,Type>& keep,
		       objective_function<Type> *of)SOURCE({
  if(!obs.hasSharedObs){
    return Type(0.0);
  }

  int nStocks = datA.size();

  vector<array<Type> > comps(datA.size());
  vector<vector<Type> > weekContrib(datA.size());
  vector<int> noYearsLAI(datA.size());
  for(int s = 0; s < datA.size(); ++s){
    array<Type> logPa = getArray(logP,s);
    comps(s) = scalePFun(confA(s), datA(s), logPa);
    weekContrib(s) = scaleWeekFun(parA(s), datA(s), logPa);
    noYearsLAI(s) = yearsPFun(confA(s),datA(s));
  }
      
  matrix<Type> predPerStock = predObsPerStock(obs, datA, confA, parA, forecastA, logF, logN, comps, logitFseason, weekContrib, mortalities, Parea, minYearAll, minAgeAll, noYearsLAI, of);
  vector<vector< MVMIX_t<Type> > > nllVecPerStock(datA.size());
  for(int s = 0; s < nllVecPerStock.size(); ++s)
    nllVecPerStock(s) = getnllVec(datA(s), confA(s), parA(s), of);

  Type nll = 0.0;

  for(int y=0;y<obs.noYears;y++){
    int totalParKey = 0;
    for(int f=0;f<obs.noFleets;f++){
      if(!((obs.fleetTypes(f)==5)||(obs.fleetTypes(f)==3)||(obs.fleetTypes(f)==6)||(obs.fleetTypes(f)==80)||(obs.fleetTypes(f)==81)||(obs.fleetTypes(f)==90))){
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
 
      }else if(obs.fleetTypes(f)==80){
	// For fleet type 80, first column of predObs gives the total prediction
	if(!isNAINT(obs.idx1(f,y))){	    
	  int iMin = obs.idx1(f,y);
	  // int iMax = dat.idx2(f,y);
	  // Setup response and prediction vectors
	  int nSeasons = CppAD::Integer(obs.auxData(iMin,5));

	  Type log_alpha = R_NegInf;
	  for(int s = 0; s < nStocks; ++s)
	    log_alpha = logspace_add_SAM(log_alpha, parA(s).logSdLogObs(confA(s).keyVarObs(f,0)));

	  vector<Type> log_X(nSeasons);
	  log_X.setConstant(R_NegInf);
	  vector<Type> log_P(nSeasons);
	  log_P.setConstant(R_NegInf);
	  // vector<Type> Keep(nSeasons);
	  // Keep.setConstant(1);
	  // vector<Type> LCDF(nSeasons);
	  // LCDF.setConstant(0);
	  // vector<Type> HCDF(nSeasons);
	  // HCDF.setConstant(0);	    
	  // for(int i=obs.idx1(f,y); i<=obs.idx2(f,y); ++i){
	  //   log_X(CppAD::Integer(obs.auxData(i,4))-1) = exp(obs.logobs(i));
	  //   log_P(CppAD::Integer(obs.auxData(i,4))-1) = predPerStock(i,0);
	  //   Keep(CppAD::Integer(obs.auxData(i,4))-1) = keep(i);
	  //   LCDF(CppAD::Integer(obs.auxData(i,4))-1) = keep.cdf_lower(i);
	  //   HCDF(CppAD::Integer(obs.auxData(i,4))-1) = keep.cdf_upper(i);
	  // }
	  data_indicator<vector<Type>, Type> K2 = keep.segment(obs.idx1(f,y),obs.idx2(f,y)-obs.idx1(f,y)+1);
	  Type xs = 1.0;
	  Type ps = 0.0;
	  for(int i=obs.idx1(f,y); i<=obs.idx2(f,y); ++i){
	    log_X(CppAD::Integer(obs.auxData(i,4))-1) = obs.logobs(i); //exp(obs.logobs(i));
	    xs += exp(obs.logobs(i));
	    log_P(CppAD::Integer(obs.auxData(i,4))-1) = predPerStock(i,0); //exp(-log_alpha + predPerStock(i,0));
	    ps += exp(predPerStock(i,0));
	  }
	  log_X(nSeasons-1) = 1.0; //log(1.0 - log_X.exp().sum());// logspace_sub_SAM(Type(0.0), logspace_sum((vector<Type>)log_X.segment(0,nSeasons-1)));
	  // log_X += (Type)0.000001;
	  // Contribution from logX -> X
	  //nll += log(xs);
	  // Contribution from X -> proportions
	  //nll -= log(fabs(obs_fun::jacobianDet(log_X)));
	  for(int i = 0; i < log_X.size(); ++i)
	    log_X(i) -= log(xs);
	  log_P(nSeasons-1) =  log(1.0 - squeeze(ps)); // exp(-log_alpha) *
	  //  vector<int> OCDF = order_keep(Keep);	      
	  // log_X(nSeasons-1) = 1.0; //log(1.0 - log_X.exp().sum());// logspace_sub_SAM(Type(0.0), logspace_sum((vector<Type>)log_X.segment(0,nSeasons-1)));	  
	  // log_X /= log_X.sum();
	  // log_X += (Type)0.000001;
	  // log_X = log_X.log();
	  // log_P(nSeasons-1) = log(1.0 - squeeze(log_P.exp().sum()));//logspace_sub_SAM(Type(0.0), logspace_sum((vector<Type>)log_P.segment(0,nSeasons-1)));
	  // From aggregation property, alpha parameters should be summed:
	  //Rcout << ddirichlet(log_X,log_P,log_alpha,Keep,true) << "\t" << ddirichlet_vtri((vector<Type>)(log_X.exp()),(vector<Type>)((log_P+log_alpha).exp()),true) << "\n";
	  nll -= ddirichlet(log_X,log_P,-log_alpha,K2,true);

	  //nll -= ddirichlet_osa(log_X,log_P,K2,true);
	  
	}
      // }else if(obs.fleetTypes(f)==81){
      // 	// For fleet type 81, first column of predObs gives the total prediction
      // 	if(!isNAINT(obs.idx1(f,y))){
      // 	  int iMin = obs.idx1(f,y);
      // 	  int iMax = obs.idx2(f,y);
      // 	  // Setup response and prediction vectors
      // 	  int nSeasons = CppAD::Integer(obs.auxData(iMin,5));
      // 	  int minAge = obs.minAgePerFleet(f);
      // 	  int maxAge = dat.maxAgePerFleet(f);
      // 	  matrix<Type> log_X(nSeasons, maxAge-minAge+1);
      // 	  log_X.setZero();
      // 	  matrix<Type> log_P(nSeasons, maxAge-minAge+1);
      // 	  log_P.setZero();
      // 	  matrix<Type> Keep(nSeasons, maxAge-minAge+1);
      // 	  Keep.setZero();
      // 	  matrix<int> nObs(nSeasons,maxAge-minAge+1);
      // 	  nObs.setZero();
      // 	  for(int i=dat.idx1(f,y); i<=dat.idx2(f,y); ++i){		
      // 	    log_X(CppAD::Integer(dat.auxData(i,4))-1, dat.aux(i,2) - minAge) = dat.logobs(i);
      // 	    log_P(CppAD::Integer(dat.auxData(i,4))-1, dat.aux(i,2) - minAge) = predPerStock(i,0);
      // 	    Keep(CppAD::Integer(dat.auxData(i,4))-1, dat.aux(i,2) - minAge) = keep(i);
      // 	    nObs(CppAD::Integer(dat.auxData(i,4))-1, dat.aux(i,2) - minAge) += 1;
      // 	  }
      // 	  // Sum to one and alpha
      // 	  vector<Type> log_alpha(maxAge-minAge+1);
      // 	  log_alpha.setConstant(R_NaReal);
      // 	  for(int i = 0; i < log_X.cols(); ++i){
      // 	    if(nObs(nSeasons-1,i)==0){
      // 	      log_X(nSeasons-1,i) = logspace_sub_SAM(Type(0.0), logspace_sum((vector<Type>)log_X.col(i).segment(0,nSeasons-1)));
      // 	      log_P(nSeasons-1,i) = logspace_sub_SAM(Type(0.0), logspace_sum((vector<Type>)log_P.col(i).segment(0,nSeasons-1)));
      // 	    }
      // 	    // From aggregation property, alpha parameters should be summed:
      // 	    log_alpha(i) = R_NegInf;
      // 	    for(int s = 0; s < nStocks; ++s)
      // 	      log_alpha(i) = parA(s).logSdLogObs(confA(s).keyVarObs(f,i + minAge - confA(s).minAge));
      // 	    nll -= ddirichlet((vector<Type>)log_X.col(i),(vector<Type>)log_P.col(i),(Type)log_alpha(i),(vector<Type>)Keep.col(i),true);
      // 	  }	      	      
      // 	}
      }else if(obs.fleetTypes(f)==90){
	if(!isNAINT(obs.idx1(f,y))){
	  int iMin = obs.idx1(f,y);
	  // int iMax = dat.idx2(f,y);
	  // Setup response and prediction vectors

	  Type log_alpha = R_NegInf;
	  for(int s = 0; s < nStocks; ++s)
	    log_alpha = logspace_add_SAM(log_alpha, parA(s).logSdLogObs(confA(s).keyVarObs(f,0)));

	  
	  vector<Type> log_X(nStocks);
	  log_X.setConstant(R_NegInf);
	  vector<Type> log_P(nStocks);
	  log_P.setConstant(R_NegInf);
	  vector<Type> Keep(nStocks);
	  Keep.setConstant(1);
	  vector<Type> LCDF(nStocks);
	  LCDF.setConstant(0);
	  vector<Type> HCDF(nStocks);
	  HCDF.setConstant(0);
	  data_indicator<vector<Type>, Type> K2 = keep.segment(obs.idx1(f,y),obs.idx2(f,y)-obs.idx1(f,y)+1);
	  Type xs = 1.0;
	  Type ps = 0.0;
	  for(int i=obs.idx1(f,y); i<=obs.idx2(f,y); ++i){
	    log_X(CppAD::Integer(obs.auxData(i,4))-1) = obs.logobs(i); //exp(obs.logobs(i));
	    xs += exp(obs.logobs(i));
	    log_P(CppAD::Integer(obs.auxData(i,4))-1) = predPerStock(i,CppAD::Integer(obs.auxData(i,4)-1)); //exp(-log_alpha + predPerStock(i,CppAD::Integer(obs.auxData(i,4)-1)));
	    ps += exp(predPerStock(i,CppAD::Integer(obs.auxData(i,4)-1)));
	    
	    Keep(CppAD::Integer(obs.auxData(i,4))-1) = keep(i);
	    LCDF(CppAD::Integer(obs.auxData(i,4))-1) = keep.cdf_lower(i);
	    HCDF(CppAD::Integer(obs.auxData(i,4))-1) = keep.cdf_upper(i);	    
	  }
	 
	  vector<int> OCDF = order_keep(Keep);
	  log_X(nStocks-1) = 0.0; //log(1.0 - log_X.exp().sum());// logspace_sub_SAM(Type(0.0), logspace_sum((vector<Type>)log_X.segment(0,nSeasons-1)));
	  // log_X += (Type)0.000001;
	  // Contribution from logX -> X
	  //nll += log(xs);
	  // Contribution from X -> proportions
	  //nll -= log(fabs(obs_fun::jacobianDet(log_X)));
	  for(int i = 0; i < log_X.size(); ++i)
	    log_X(i) -= log(xs);
	  log_P(nStocks-1) = log((1.0 - squeeze(ps))); //exp(-log_alpha) * (1.0 - squeeze(ps));	     
	  // log_X /= log_X.sum();
	  // log_X = log_X.log();
	  // log_P(nStocks-1) = log(1.0 - squeeze(log_P.exp().sum())); //logspace_sub_SAM(Type(0.0), logspace_sum((vector<Type>)log_P.segment(0,nStocks-1)));
	  // Rcout << log_X.exp().sum() << "\n";
	  // Rcout << log_P.exp().sum() << "\n";
	  // Rcout << log_alpha << "\n";
	  // Rcout << ddirichlet(log_X,log_P,log_alpha,Keep,true) << "\t" << ddirichlet_vtri((vector<Type>)(log_X.exp()),(vector<Type>)((log_P+log_alpha).exp()),true) << "\n\n";
	  // Rcout << "From 90\n";
	  // Rcout << "log_X:\n" << log_X << "\n";
	  // Rcout << "log_P:\n" << log_P << "\n";
	  // Rcout << "log_alpha:\n" << log_alpha << "\n";
	  // Rcout << "Keep:\n" << Keep << "\n";
	  // Rcout << ddirichlet(log_X,log_P,0*log_alpha,Keep,LCDF,HCDF,OCDF,true) << "\n";
	  // Rcout << ddirichlet_osa(log_X,log_P,K2,true) << "\n";
	  nll -= ddirichlet(log_X,log_P,-log_alpha,K2,true);
	  //nll -= ddirichlet_osa(log_X,log_P,K2,true);
	  // for(int j = 0; j < log_X.size(); ++j)
	  //   nll += ((log_X(j) - log_P(j)) * (log_X(j) - log_P(j)));
	  //nll -= ddirichlet_vtri((vector<Type>)(log_X.exp()),(vector<Type>)((log_P+log_alpha).exp()),true);
	}
	     }else if(obs.fleetTypes(f)==90){
	if(!isNAINT(obs.idx1(f,y))){
	  int iMin = obs.idx1(f,y);
	  // int iMax = dat.idx2(f,y);
	  // Setup response and prediction vectors
	  Type log_alpha = R_NegInf;
	  for(int s = 0; s < nStocks; ++s)
	    log_alpha = logspace_add_SAM(log_alpha, parA(s).logSdLogObs(confA(s).keyVarObs(f,0)));

	  
	  vector<Type> log_X(nStocks);
	  log_X.setConstant(R_NegInf);
	  vector<Type> log_P(nStocks);
	  log_P.setConstant(R_NegInf);
	  data_indicator<vector<Type>, Type> K2 = keep.segment(obs.idx1(f,y),obs.idx2(f,y)-obs.idx1(f,y)+1);
	  Type xs = 1.0;
	  Type ps = 0.0;
	  for(int i=obs.idx1(f,y); i<=obs.idx2(f,y); ++i){
	    log_X(CppAD::Integer(obs.auxData(i,4))-1) = obs.logobs(i); //exp(obs.logobs(i));
	    xs += exp(obs.logobs(i));
	    log_P(CppAD::Integer(obs.auxData(i,4))-1) = predPerStock(i,0); //exp(-log_alpha + predPerStock(i,CppAD::Integer(obs.auxData(i,4)-1)));
	    ps += exp(predPerStock(i,0));	    
	  }
	 
	  log_X(nStocks-1) = 0.0;
	  for(int i = 0; i < log_X.size(); ++i)
	    log_X(i) -= log(xs);
	  log_P(nStocks-1) = log((1.0 - squeeze(ps))); //exp(-log_alpha) * (1.0 - squeeze(ps));	     
	  nll -= ddirichlet(log_X,log_P,-log_alpha,K2,true);	  
	}
      }else{
	Rf_error("fleetType not implemented for shared observations");
      }
    }
  }
  return nll; //nll;
  
			 });


MSM_SPECIALIZATION(double sharedObservation(shared_obs<double>&, vector<dataSet<double> >&, vector<confSet>&, vector<paraSet<double> >&, vector<forecastSet<double> >&, cmoe_matrix<double>&, cmoe_matrix<double>&, cmoe_matrix<double>&, cmoe_3darray<double>&, vector<MortalitySet<double> >&, matrix<double>&, int, int, data_indicator<vector<double>, double>&, objective_function<double>*));
MSM_SPECIALIZATION(TMBad::ad_aug sharedObservation(shared_obs<TMBad::ad_aug>&, vector<dataSet<TMBad::ad_aug> >&, vector<confSet>&, vector<paraSet<TMBad::ad_aug> >&, vector<forecastSet<TMBad::ad_aug> >&, cmoe_matrix<TMBad::ad_aug>&, cmoe_matrix<TMBad::ad_aug>&, cmoe_matrix<TMBad::ad_aug>&, cmoe_3darray<TMBad::ad_aug>&, vector<MortalitySet<TMBad::ad_aug> >&, matrix<TMBad::ad_aug>&, int, int, data_indicator<vector<TMBad::ad_aug>, TMBad::ad_aug>&, objective_function<TMBad::ad_aug>*));
			 
