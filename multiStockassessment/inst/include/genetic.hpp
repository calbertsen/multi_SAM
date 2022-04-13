
template<class Type>
vector<Type> toLogProportionG(vector<Type> x){
  vector<Type> r(x.size()+1);
  r.setZero();
  Type lps = 0.0;
  for(int i = 0; i < x.size(); ++i){
    r(i) = x(i);
    lps = logspace_add2(lps,x(i));
  }
  r -= lps;
  return r;
}


template<class Type>
Type dmultinom2(vector<Type> x, vector<Type> logp, int give_log=0){
  Type logres = 0.0;
  Type xsum = 0.0;
  for(int i = 0; i < x.size(); ++i){
    logres += -lgamma(x(i) + 1.0) + x(i) * logp(i);
    xsum += x(i);
  }
  logres += lgamma(xsum + 1.0);
  //-lgamma(alpha.sum()) + lgamma(alpha).sum() + xcontr.sum();
  if(give_log)return logres; else return exp(logres);
}


template<class Type>
Type ddirichletmultinom(vector<Type> x, vector<Type> logp, Type logAlphaScale, int give_log=0){
  vector<Type> alpha = exp(logp + logAlphaScale);
  Type logres = 0.0;
  Type alphsum = 0.0;
  Type xsum = 0.0;
  for(int i = 0; i < x.size(); ++i){
    logres += lgamma(x(i) + alpha(i)) - lgamma(x(i) + 1.0) - lgamma(alpha(i));
    alphsum += alpha(i);
    xsum += x(i);
  }
  logres += lgamma(xsum + 1.0) + lgamma(alphsum) - lgamma(xsum + alphsum);
  //-lgamma(alpha.sum()) + lgamma(alpha).sum() + xcontr.sum();
  if(give_log)return logres; else return exp(logres);
}



// x:  #Allele x # Loci matrix of observations
// mu:  #(Allele-1) x #Loci matrix of allele frequencies ; 
// scale: vector of length #Loci; 
template<class Type>
Type dAlleleCount(matrix<Type> x, matrix<Type> mu, vector<Type> scale, int give_log = 0){
  Type r = 0.0;
  for(int i = 0; i < x.cols(); ++i){
    vector<Type> logp = toLogProportionG((vector<Type>)mu.col(i));
    vector<Type> y = x.col(i);
    if(sum(y) == 0){
      r += 0.0;
    }else if(R_finite(asDouble(scale(i)))){
      r += ddirichletmultinom(y,logp,scale(i),true);
    }else if(!R_finite(asDouble(scale(i)))){
    }else{
      r += dmultinom2(y,logp,true);
    }
  }
  if(give_log)return r; else return exp(r);
}



template<class Type>
struct genetic_sample {
  int keyStock;			// NA if not baseline
  int keyGridTime;
  int keyGridSpace;
  int keyTrip;
  int year;
  int age;
  int fleet;
  matrix<Type> alleleCount;	// Matrix of dim nAllele x nLoci

  // Constructors
  genetic_sample(){};
  genetic_sample(SEXP x){
    keyStock = (int)*REAL(getListElement(x,"keyStock"));
    keyGridTime = (int)*REAL(getListElement(x,"keyGridTime"));
    keyGridSpace = (int)*REAL(getListElement(x,"keyGridSpace"));
    keyTrip = (int)*REAL(getListElement(x,"keyTrip"));
    year = (int)*REAL(getListElement(x,"year"));
    age = (int)*REAL(getListElement(x,"age"));
    fleet = (int)*REAL(getListElement(x,"fleet"));
    alleleCount = asMatrix<Type>(getListElement(x,"alleleCount"));      
  }
};


template<class Type>
struct genetic_data {
  vector<genetic_sample<Type> > samples;
  Eigen::SparseMatrix<Type> Qspace;
  Eigen::SparseMatrix<Type> Qtime;
  int QorderTime;
  int QorderSpace;
  matrix<Type> stock2gen;	// nStocks x nGenetics
  // Constructors
  genetic_data() : samples(0) {};
  genetic_data(SEXP x){
    SEXP samp = PROTECT(getListElement(x,"samples"));
    int n = Rf_length(samp);
    samples = vector<genetic_sample<Type> >(n);
    for(int i = 0; i < n; ++i){
      SEXP y = PROTECT(VECTOR_ELT(samp,i));
      samples(i) = genetic_sample<Type>(y);
      UNPROTECT(1);
    }
    UNPROTECT(1);
    Qspace = tmbutils::asSparseMatrix<Type>(getListElement(x,"Qspace"));
    Qtime = tmbutils::asSparseMatrix<Type>(getListElement(x,"Qtime"));
    QorderTime = (int)*REAL(getListElement(x,"QorderTime"));
    QorderSpace = (int)*REAL(getListElement(x,"QorderSpace"));
    stock2gen = asMatrix<Type>(getListElement(x,"stock2gen"));
  }
};


template<class Type>
struct genetic_parameters {
  Type logKappaSpace;
  Type logKappaTime;
  vector<Type> corparST;
  vector<Type> logSdST;
  Type corparAge;
  vector<Type> corparTrip;
  vector<Type> logSdTrip;
  array<Type> alleleFreq;// (nAllele-1) x nLoci x nStockGenetic
  matrix<Type> dmScale;
  matrix<Type> muLogP;		// nAge x nStockGenetic
};


template<class Type>
Type nllGenetics(shared_obs<Type>& obs,
		 vector<dataSet<Type> >& datA,
		 vector<confSet>& confA,
		 vector<paraSet<Type> >& parA,
		 vector<forecastSet<Type> >& forecastA,
		 cmoe_matrix<Type>& logF,
		 cmoe_matrix<Type>& logN,
		 genetic_parameters<Type>& genpar,
		 genetic_data<Type>& gendat,
		 //array<Type> alleleFreq, // (nAllele-1) x nLoci x nStockGenetic
		 array<Type>& logGst, // nSpace x nTime x nAges x (nStock-1) (??order)
		 matrix<Type>& logGtrip, // (nStock-1) x nTrips
		 int maxAgeAll,
		 int minAgeAll
		 ){

  if(gendat.samples.size() == 0)
    return 0.0;
  int nStockManagement = gendat.stock2gen.rows();
  int nStockGenetic = gendat.stock2gen.cols();
  // int nre_Space = logGst.dim[0];
  // int nre_Time = logGst.dim[1];
  int nre_Age = logGst.dim[2];
  // int nre_Stock = logGst.dim[3];
  Type nll = 0.0;
  // nll form space-time-stock random effects
  using namespace density;
  Eigen::SparseMatrix<Type> Qspace = gendat.Qspace;
  for(int i = 0; i < Qspace.rows(); ++i)
    Qspace.coeffRef(i,i) += exp(genpar.logKappaSpace);
  Eigen::SparseMatrix<Type> Qtime = gendat.Qtime;
  for(int i = 0; i < Qtime.rows(); ++i)
    Qtime.coeffRef(i,i) += exp(genpar.logKappaTime);

  if(Qtime.cols() > 0 && Qspace.cols() > 0)
    nll += SEPARABLE(SEPARABLE(SEPARABLE(
					 VECSCALE(UNSTRUCTURED_CORR(genpar.corparST),exp(genpar.logSdST)), // Stock
					 AR1(invlogit(genpar.corparAge))), // Age
			       GMRF(Qtime, gendat.QorderTime)), // Time
		     GMRF(Qspace, gendat.QorderSpace))(logGst); // Space
  // nll for trip random effects
  VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > nllTrip = VECSCALE(UNSTRUCTURED_CORR(genpar.corparTrip),exp(genpar.logSdTrip));
  for(int i = 0; i < logGtrip.cols(); ++i)
    nll += nllTrip((vector<Type>)logGtrip.col(i));
  // nll for genetic samples
  for(int i = 0; i < gendat.samples.size(); ++i){
    genetic_sample<Type> gs = gendat.samples(i);
    // Type llTmp = R_NegInf;
    // Stock probabilities
    vector<Type> logPnn(nStockGenetic);
    logPnn.setConstant(R_NegInf);
    if(gs.keyStock == R_NaInt){	// Stock is unknown (i.e. new informative sample)
        if(obs.hasSharedObs == 0)
	  Rf_error("Informative genetic samples require shared observations");

      vector<Type> predCatch(nStockGenetic); // nGeneticStocks
      predCatch.setConstant(R_NegInf);
      int a1 = minAgeAll;
      int a2 = maxAgeAll;
      if(gs.age != R_NaInt){	// if age is known, reduce loop
	a1 = gs.age;
	a2 = gs.age;
      }

      for(int aa = a1; aa <= a2; ++aa){
	vector<Type> PCtmpM(nStockManagement); // Should be nModelledStocks (not nGeneticStocks)
	PCtmpM.setZero();
	vector<Type> PCtmpG(nStockGenetic); // Should be nModelledStocks (not nGeneticStocks)
	PCtmpG.setZero();
	if(gs.fleet == R_NaInt){ // No fleet, corresponds to total stock at time of obs
	  Rf_error("Missing fleet is not implemented yet for genetic data");
	}else{
	  PCtmpM = predOneObsPerStock(gs.fleet,
				     obs.fleetTypes(gs.fleet-1),
				     aa,
				     gs.year,
				     datA,
				     confA,
				     parA,
				      forecastA,
				     logF,
				     logN,
				     (vector<Type>)obs.keyFleetStock.row(gs.fleet-1));
	}
	// Map modelled stocks to genetic stocks
	PCtmpG = PCtmpM.matrix().transpose() * gendat.stock2gen;
	for(int s = 0; s < nStockGenetic; ++s){
	  Type PCcorrected = PCtmpG(s) + genpar.muLogP(aa - minAgeAll,s);
	  if(s < nStockGenetic-1){
	    if(nre_Age == maxAgeAll - minAgeAll + 1){
	      if(gs.keyGridSpace != R_NaInt && gs.keyGridTime != R_NaInt)
		PCcorrected += logGst(gs.keyGridSpace, gs.keyGridTime, aa - minAgeAll, s);
	    }else{
	      if(gs.keyGridSpace != R_NaInt && gs.keyGridTime != R_NaInt)
		PCcorrected += logGst(gs.keyGridSpace, gs.keyGridTime, 0, s);
	    }
	    if(gs.keyTrip != R_NaInt){
	      PCcorrected += logGtrip(s, gs.keyTrip);
	    }
	  }
	  predCatch(s) = logspace_add2(predCatch(s), PCcorrected);
	}	
      }
      // if(nre_Age == 1 && (maxAgeAll - minAgeAll > 0)){
      // 	for(int s = 0; s < nStockGenetic-1; ++s)
      // 	  predCatch(s) += logGst(gs.keyGridSpace, gs.keyGridTime, 0, s);
      // }
      // if(gs.keyTrip != R_NaInt)
      // 	for(int s = 0; s < nStock-1; ++s)
      // 	  predCatch(s) += logGtrip(s, gs.keyTrip);
      logPnn = predCatch;
    }else{// Stock is known (i.e. baseline)
      logPnn(gs.keyStock) = 0.0;
    }
    Type logPZ = logspace_sum(logPnn);
    vector<Type> logP = logPnn - logPZ;

    // Data contribution
    Type llTmp = R_NegInf;
    for(int s = 0; s < nStockGenetic; ++s){
      Type n0 = dAlleleCount(gs.alleleCount,
			     genpar.alleleFreq.col(s).matrix(),
			     (vector<Type>)genpar.dmScale.col(s),
			     true);
      llTmp = logspace_add2(llTmp, n0 + logP(s));      
    }
    nll -= llTmp;
  }
  return nll;
}
