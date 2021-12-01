

template<class Type>
vector<Type> predOneObsPerStock(int fleet,	// obs.aux(i,1)
				int fleetType,	// obs.fleetTypes(f-1)
				int age,	// obs.aux(i,2)-confA(s).minAge
				int year, // obs.aux(i,0)
				vector<dataSet<Type> > datA,
				vector<confSet> confA,
				vector<paraSet<Type>> parA,
				cmoe_matrix<Type> logF,
				cmoe_matrix<Type> logN,
				vector<Type> keyFleetStock){
  if(fleetType > 2)
    Rf_error("fleetType not implemented for shared observations");
  vector<Type> pred(datA.size());
  pred.setConstant(R_NegInf);
  for(int s = 0; s < datA.size(); ++s){
    // Skip if age or year is not used for this stock
    int y = year -  CppAD::Integer(datA(s).years(0));
    int a = age - confA(s).minAge;
    if(y < 0 || y >= datA(s).noYears + datA(s).forecast.nYears)
      continue;
    if(a < 0 || a > (confA(s).maxAge - confA(s).minAge))
      continue;

    array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
    logNa = logN.col(s);
    array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
    logFa = logF.col(s);    
    pred(s) = predOneObs(fleet,	// obs.aux(i,1)
			 fleetType,	// obs.fleetTypes(f-1)
			 age,	// obs.aux(i,2)-confA(s).minAge
			 year, // obs.aux(i,0)
			 CppAD::Integer(datA(s).years(0)), // minYear
			 datA(s),
			 confA(s),
			 parA(s),
			 logFa,
			 logNa,
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
}

template<class Type>
matrix<Type> predObsPerStock(shared_obs<Type> obs,
		       vector<dataSet<Type> > datA,
		       vector<confSet> confA,
		       vector<paraSet<Type> > parA,
		       cmoe_matrix<Type> logF,
		       cmoe_matrix<Type> logN,
		       // vector<Type> logSdObs,
		       int minYearAll,
		       int minAgeAll,
		       objective_function<Type> *of){

  matrix<Type> predPerStock(obs.aux.dim[0], datA.size());
  predPerStock.setConstant(R_NegInf);

  // for(int s = 0; s < datA.size(); ++s){
  //   int f, ft, a, y, ma, pg;
  //   //, yy, scaleIdx, ma, pg;  // a is no longer just ages, but an attribute (e.g. age or length) 
  //   Type logzz=Type(R_NegInf);
    for(int i=0;i<obs.aux.dim[0];i++){
      // array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
      // logNa = logN.col(s);
      // array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
      // logFa = logF.col(s);

      // y=obs.aux(i,0) - CppAD::Integer(datA(s).years(0));
      // f=obs.aux(i,1);
      // ft=obs.fleetTypes(f-1);
      // a=obs.aux(i,2)-confA(s).minAge;

      predPerStock.row(i) = predOneObsPerStock(obs.aux(i,1), // fleet
					       obs.fleetTypes(obs.aux(i,1)-1), // fleetType
					       obs.aux(i,2), // age
					       obs.aux(i,0), // year
					       datA,
					       confA,
					       parA,
					       logF,
					       logN,
					       (vector<Type>)obs.keyFleetStock.row(obs.aux(i,1)-1));
    }
					       
  //     // Skip if age or year is not used for this stock
  //     if(y < 0 || y >= datA(s).noYears + datA(s).forecast.nYears)
  // 	continue;
  //     if(a < 0 || a > (confA(s).maxAge - confA(s).minAge))
  // 	continue;

  //     if(obs.aux(i,2)==datA(s).maxAgePerFleet(f-1)){ma=1;}else{ma=0;}
  //     pg=confA(s).maxAgePlusGroup(f-1);
      
  //     // Set log(Z)
  //     if(ft<3){ 
  // 	logzz = log(datA(s).natMor(y,a));
  // 	if(confA(s).keyLogFsta(0,a)>(-1)){
  // 	  logzz = logspace_add2(logzz, logFa(confA(s).keyLogFsta(0,a),y));
  // 	}
  //     }    

  //     // Predict based on fleet type
  //     switch(ft){
  //     case 0:
  // 	//pred(i)=logN(a,y)-logzz+log(1-exp(-exp(logzz)));
  // 	predPerStock(i,s)=logNa(a,y)-logzz+logspace_sub2(Type(0.0),-exp(logzz));
  // 	if(confA(s).keyLogFsta(f-1,a)>(-1)){
  // 	  predPerStock(i,s)+=logFa(confA(s).keyLogFsta(0,a),y);
  // 	}
  // 	// scaleIdx=-1;
  // 	// yy=dat.aux(i,0);
  // 	// for(int j=0; j<conf.noScaledYears; ++j){
  // 	//   if(yy==conf.keyScaledYears(j)){
  // 	//     scaleIdx=conf.keyParScaledYA(j,a);
  // 	//     if(scaleIdx>=0){
  // 	//       pred(i)-=par.logScale(scaleIdx);
  // 	//     }
  // 	//     break;
  // 	//   }
  // 	// }
  // 	break;
  //     case 2:
  // 	if((pg!=confA(s).maxAgePlusGroup(0))&&(a==(confA(s).maxAge-confA(s).minAge))){
  //         Rf_error("When maximum age for the fleet is the same as maximum age in the assessment it must be treated the same way as catches w.r.t. plusgroup configuration");
  // 	}

  // 	if((ma==1) && (pg==1)){
  // 	  for(int aa=a; aa<=(confA(s).maxAge-confA(s).minAge); aa++){
  // 	    logzz = log(datA(s).natMor(y,aa));
  //           if(confA(s).keyLogFsta(0,aa)>(-1)){
  //             logzz = logspace_add2(logzz, logFa(confA(s).keyLogFsta(0,aa),y));
  // 	    }
  // 	    predPerStock(i,a) = logspace_add2(predPerStock(i,a), logNa(aa,y)-exp(logzz)*datA(s).sampleTimes(f-1));
  // 	  }
  // 	}else{
  //         predPerStock(i,s)=logNa(a,y)-exp(logzz)*datA(s).sampleTimes(f-1);
  // 	}
  //       if(confA(s).keyQpow(f-1,a)>(-1)){
  //         predPerStock(i,s)*=exp(parA(s).logQpow(confA(s).keyQpow(f-1,a))); 
  //       }
  //       if(confA(s).keyLogFpar(f-1,a)>(-1)){
  //         predPerStock(i,s)+=parA(s).logFpar(confA(s).keyLogFpar(f-1,a));
  //       }
  // 	break;
  //     default:
  // 	Rf_error("Unknown fleet code");
  // 	break;
  //     }
  //     if(obs.keyFleetStock(f-1,s) == 0){
  // 	predPerStock(i,s) = R_NegInf;
  // 	//logSdPerStock(i,s) = R_NegInf;      
  //     }else{
  // 	predPerStock(i,s) += log(obs.keyFleetStock(f-1,s));
  // 	//logSdPerStock(i,s) += log(obs.keyFleetStock(f-1,s));
  //     }
  //   } // End loop over observtions
  // }   // End loop over stocks
  // Sum catch predictions
  // vector<Type> predObsShared(predPerStock.rows());
  // predObsShared.setConstant(R_NegInf);
  // vector<Type> varObsShared(logSdPerStock.rows());
  // varObsShared.setConstant(R_NegInf);

  // for(int i = 0; i < predPerStock.rows(); ++i){
  //   vector<Type> loggg(predPerStock.cols());
  //   loggg.setZero();
  //   for(int s = 0; s < predPerStock.cols(); ++s){
  //     predObsShared(i) = logspace_add2(predObsShared(i), predPerStock(i,s));
  //     loggg(s) = predPerStock(i,s);
  //   }
  //   loggg -= predObsShared(i);
  //   for(int s = 0; s < predPerStock.cols(); ++s){
  //     varObsShared(i) = logspace_add2(varObsShared(i), 2.0 * loggg(s) + logSdPerStock(i,s));
  //   }    
  // }

  return predPerStock;
}




template<class Type>
Type sharedObservation(shared_obs<Type> obs,
		       vector<dataSet<Type> > datA,
		       vector<confSet> confA,
		       vector<paraSet<Type> > parA,
		       cmoe_matrix<Type> logF,
		       cmoe_matrix<Type> logN,
		       // vector<Type> logSdObs,
		       int minYearAll,
		       int minAgeAll,
		       data_indicator<vector<Type>,Type> keep,
		       objective_function<Type> *of){

  if(!obs.hasSharedObs){
    return Type(0.0);
  }
  matrix<Type> predPerStock = predObsPerStock(obs, datA, confA, parA, logF, logN, minYearAll, minAgeAll, of);
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
	      predObsSegment(a) = logspace_add2(predObsSegment(a), predUse(a,s));
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
              
	    // for(int idxV=0; idxV<currentVar.size(); ++idxV){
	    //   if(isNA(dat.weight(idxfrom+idxV))){
	    // 	sqrtW(idxV)=Type(1.0);
	    // 	int a = dat.aux(idxfrom+idxV,2)-conf.minAge;
	    // 	if(conf.predVarObsLink(f,a)>(-1)){
	    // 	  sqrtW(idxV) = sqrt(findLinkV(par.logSdLogObs(conf.keyVarObs(f,a))+(exp(par.predVarObs(conf.predVarObsLink(f,a))) -Type(1))*predObs(idxfrom+idxV))/currentVar(idxV));
	    // 	}
	    //   }else{
	    // 	if(conf.fixVarToWeight==1){
	    // 	  sqrtW(idxV)=sqrt(dat.weight(idxfrom+idxV)/currentVar(idxV));
	    // 	}else{
	    // 	  sqrtW(idxV)=sqrt(Type(1)/dat.weight(idxfrom+idxV));
	    // 	}
	    //   }
	    // }
	    // if(isNAINT(dat.idxCor(f,y))){
	    nll += thisNll((obs.logobs.segment(idxfrom,idxlength)-predObsSegment), keep.segment(idxfrom,idxlength)); ///sqrtW,keep.segment(idxfrom,idxlength));
	    //     nll += (log(sqrtW)*keep.segment(idxfrom,idxlength)).sum();
	      SIMULATE_F(of){
		obs.logobs.segment(idxfrom,idxlength) = predObsSegment + (thisNll.simulate()*sqrtW);
	      }
	    // }else{
	    //   int thisdim=currentVar.size();
	    //   matrix<Type> thiscor=dat.corList(dat.idxCor(f,y));
	    //   matrix<Type> thiscov(thisdim,thisdim);
	    //   for(int r=0;r<thisdim;++r){
	    // 	for(int c=0;c<thisdim;++c){
	    // 	  thiscov(r,c)=thiscor(r,c)*sqrt(currentVar(r)*currentVar(c));
	    // 	}
	    //   } 
	    //   MVMIX_t<Type> thisnll(thiscov,conf.fracMixObs(f));
	    //   nll+= thisnll((obs.logobs.segment(idxfrom,idxlength)-predObsSegment)/sqrtW, keep.segment(idxfrom,idxlength));              
	    //   nll+= (log(sqrtW)*keep.segment(idxfrom,idxlength)).sum();
	    //   SIMULATE_F(of){
	    //     dat.logobs.segment(idxfrom,idxlength) = predObs.segment(idxfrom,idxlength) + thisnll.simulate()*sqrtW;
	    //   }
            // }
	    break;
	  case 1: // (ALN) Additive logistic-normal proportions + log-normal total numbers
	    nll +=  thisNll(addLogratio((vector<Type>)obs.logobs.segment(idxfrom,idxlength))-addLogratio((vector<Type>)predObsSegment));
	    nll += log(log2proportion((vector<Type>)obs.logobs.segment(idxfrom,idxlength))).sum();
	    nll -= dnorm(log(log2expsum((vector<Type>)obs.logobs.segment(idxfrom,idxlength))),
	     	         log(log2expsum((vector<Type>)predObsSegment)),
	   	         exp(parA(0).logSdLogTotalObs(totalParKey++)),true); // SHOULD BE CHANGED TO USE DELTA METHOD
	    nll += log(log2expsum((vector<Type>)obs.logobs.segment(idxfrom,idxlength)));
	    nll -= log(fabs(jacobianDet((vector<Type>)obs.logobs.segment(idxfrom,idxlength).exp())));
            nll -= obs.logobs.segment(idxfrom,idxlength).sum();
	    SIMULATE_F(of){
	      vector<Type> logProb(idxlength);
	      logProb.setZero();
	      logProb.segment(0,idxlength-1) = addLogratio(((vector<Type>)predObsSegment)) + thisNll.simulate();
	      Type logDenom = logExpSum(logProb);
	      logProb -= logDenom;
	      Type logTotal = rnorm(log(log2expsum((vector<Type>)predObsSegment)),
	  			    exp(parA(0).logSdLogTotalObs(totalParKey++)));
	      obs.logobs.segment(idxfrom,idxlength) = logProb + logTotal; 
	    }
	    break;
	  default:
	    Rf_error("Unknown obsLikelihoodFlag");
	  }
        }

      // }else{ //dat.fleetTypes(f)==5
      //   if(dat.fleetTypes(f)==5){
      //     if(!isNAINT(dat.idx1(f,y))){    
      //       for(int i=dat.idx1(f,y); i<=dat.idx2(f,y); ++i){
      //         nll += -keep(i)*dnbinom(dat.logobs(i),predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i),true);
      //         SIMULATE_F(of){
      // 		dat.logobs(i) = rnbinom(predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i));
      //         }
      //       }
      //     }
      //   }else{
      //     if(dat.fleetTypes(f)==3){
      //       Type sd=0;
      //       if(!isNAINT(dat.idx1(f,y))){
      //         for(int i=dat.idx1(f,y); i<=dat.idx2(f,y); ++i){
      //           if(conf.keyBiomassTreat(f)==3){
      //             sd = sqrt(varLogCatch(y));
      //           }else{
      //             if(conf.keyBiomassTreat(f)==4){
      //               sd = sqrt(varLogLand(y));
      //             }else{
      //               sd = exp(par.logSdLogObs(conf.keyVarObs(f,0)));
      //             }
      //           }  
      //           nll += -keep(i)*dnorm(dat.logobs(i),predObs(i),sd,true);
      //           SIMULATE_F(of){
      // 		  dat.logobs(i) = rnorm(predObs(i),sd);
      //           }
      //         }
      //       }    
      //     }
        // }   
      }   
    }  
  }
  return nll; //nll;
}


