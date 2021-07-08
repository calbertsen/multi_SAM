

template<class Type>
Type sharedObservation(shared_obs<Type> obs,
		       vector<dataSet<Type> > datA,
		       vector<confSet> confA,
		       vector<paraSet<Type> > parA,
		       cmoe_matrix<Type> logF,
		       cmoe_matrix<Type> logN,
		       vector<Type> logSdObs,
		       int minYearAll,
		       int minAgeAll,
		       objective_function<Type> *of){

  if(!obs.hasSharedObs){
    return Type(0.0);
  }
  
 
  // Generate catch predictions for each stock involved
  // logF/logN could/should be scaled if the combined catch only covers part of the stock
  matrix<Type> predPerStock(obs.aux.dim[0], datA.size());
  predPerStock.setConstant(R_NegInf);

  matrix<Type> logSdPerStock(obs.aux.dim[0], datA.size());
  logSdPerStock.setConstant(R_NegInf);

  

  for(int s = 0; s < datA.size(); ++s){
    int f, ft, a, y, ma, pg;
    //, yy, scaleIdx, ma, pg;  // a is no longer just ages, but an attribute (e.g. age or length) 
    Type logzz=Type(R_NegInf);
    for(int i=0;i<obs.aux.dim[0];i++){
      array<Type> logNa(logN.col(s).rows(),logN.col(s).cols());
      logNa = logN.col(s);
      array<Type> logFa(logF.col(s).rows(),logF.col(s).cols());
      logFa = logF.col(s);

      y=obs.aux(i,0) - CppAD::Integer(datA(s).years(0));
      f=obs.aux(i,1);
      ft=obs.fleetTypes(f-1);
      a=obs.aux(i,2)-confA(s).minAge;
      
      // Skip if age or year is not used for this stock
      if(y < 0 || y >= datA(s).noYears + datA(s).forecast.nYears)
	continue;
      if(a < 0 || a > (confA(s).maxAge - confA(s).minAge))
	continue;

      if(obs.aux(i,2)==datA(s).maxAgePerFleet(f-1)){ma=1;}else{ma=0;}
      pg=confA(s).maxAgePlusGroup(f-1);
      
      // Set log(Z)
      if(ft<3){ 
	logzz = log(datA(s).natMor(y,a));
	if(confA(s).keyLogFsta(0,a)>(-1)){
	  logzz = logspace_add2(logzz, logFa(confA(s).keyLogFsta(0,a),y));
	}
      }    

      // Variance
      logSdPerStock(i,s) = parA(s).logSdLogObs(confA(s).keyVarObs(f-1,a));
      // Predict based on fleet type
      switch(ft){
      case 0:
	//pred(i)=logN(a,y)-logzz+log(1-exp(-exp(logzz)));
	predPerStock(i,s)=logNa(a,y)-logzz+logspace_sub2(Type(0.0),-exp(logzz));
	if(confA(s).keyLogFsta(f-1,a)>(-1)){
	  predPerStock(i,s)+=logFa(confA(s).keyLogFsta(0,a),y);
	}
	// scaleIdx=-1;
	// yy=dat.aux(i,0);
	// for(int j=0; j<conf.noScaledYears; ++j){
	//   if(yy==conf.keyScaledYears(j)){
	//     scaleIdx=conf.keyParScaledYA(j,a);
	//     if(scaleIdx>=0){
	//       pred(i)-=par.logScale(scaleIdx);
	//     }
	//     break;
	//   }
	// }
	break;
      case 2:
	if((pg!=confA(s).maxAgePlusGroup(0))&&(a==(confA(s).maxAge-confA(s).minAge))){
          Rf_error("When maximum age for the fleet is the same as maximum age in the assessment it must be treated the same way as catches w.r.t. plusgroup configuration");
  	}

	if((ma==1) && (pg==1)){
	  for(int aa=a; aa<=(confA(s).maxAge-confA(s).minAge); aa++){
	    logzz = log(datA(s).natMor(y,aa));
            if(confA(s).keyLogFsta(0,aa)>(-1)){
              logzz = logspace_add2(logzz, logFa(confA(s).keyLogFsta(0,aa),y));
	    }
	    predPerStock(i,a) = logspace_add2(predPerStock(i,a), logNa(aa,y)-exp(logzz)*datA(s).sampleTimes(f-1));
	  }
	}else{
          predPerStock(i,s)=logNa(a,y)-exp(logzz)*datA(s).sampleTimes(f-1);
	}
        if(confA(s).keyQpow(f-1,a)>(-1)){
          predPerStock(i,s)*=exp(parA(s).logQpow(confA(s).keyQpow(f-1,a))); 
        }
        if(confA(s).keyLogFpar(f-1,a)>(-1)){
          predPerStock(i,s)+=parA(s).logFpar(confA(s).keyLogFpar(f-1,a));
        }
	break;
      default:
	Rf_error("Unknown fleet code");
	return 0 ;
	break;
      }
      if(obs.keyFleetStock(f-1,s) == 0){
	predPerStock(i,s) = R_NegInf;
	logSdPerStock(i,s) = R_NegInf;      
      }else{
	predPerStock(i,s) += log(obs.keyFleetStock(f-1,s));
	logSdPerStock(i,s) += log(obs.keyFleetStock(f-1,s));
      }
    } // End loop over observtions
  }   // End loop over stocks
  // Sum catch predictions
  vector<Type> predObsShared(predPerStock.rows());
  predObsShared.setConstant(R_NegInf);
  vector<Type> logVarObsShared(logSdPerStock.rows());
  logVarObsShared.setConstant(R_NegInf);

  for(int i = 0; i < predPerStock.rows(); ++i)
    for(int j = 0; j < predPerStock.cols(); ++j){
      predObsShared(i) = logspace_add2(predObsShared(i), predPerStock(i,j));
      logVarObsShared(i) = logspace_add2(logVarObsShared(i), 2.0 * logSdPerStock(i,j));
    }

  REPORT_F(predPerStock,of);
  REPORT_F(predObsShared,of);
  // Prepare covariance

  // Return likelihood contribution
  Type nll = 0.0;
  for(int i = 0; i < obs.logobs.size(); ++i){
    //int f=obs.aux(i,1);
    nll -= dnorm((Type)obs.logobs(i), (Type)predObsShared(i), (Type)exp(0.5 * logVarObsShared(i)), true);
  }

  return nll; //nll;
}
