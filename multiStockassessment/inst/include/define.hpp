
// Define data structure

template<class Type>
struct sam_data {

  vector<dataSet<Type> > dataSets;
  vector<confSet> confSets;

  sam_data(SEXP x) : dataSets(Rf_length(x)), confSets(Rf_length(x)){ // Constructor
    int n = Rf_length(x);

    // Loop through x to get the data for each area
    for(int i = 0; i < n; ++i){
      SEXP y = PROTECT(VECTOR_ELT(x,i));
      dataSets(i) = dataSet<Type>(y);
      confSets(i) = confSet(y);
      UNPROTECT(1);
    }
  }; // End constructor

};




template<class Type>
struct shared_obs {

  enum CovCombineType {
		       CC_Delta_LogSumExp,
		       CC_Delta_SumExp,
		       CC_Average
  };
  
  int hasSharedObs;
  vector<int> fleetTypes;
  vector<int> maxAgePlusGroup;
  array<int> aux;
  vector<Type> logobs;
  matrix<Type> keyFleetStock;
  int noYears;
  int noFleets;
  array<int> idx1;
  array<int> idx2;
  array<int> idxCor;

  vector<CovCombineType> covCombine;

  shared_obs() : hasSharedObs(0),
		 fleetTypes(),
		 maxAgePlusGroup(),
		 aux(),
		 logobs(),
		 keyFleetStock(),
		 noYears(0),
		 noFleets(0),
		 idx1(),
		 idx2(),
		 idxCor(),
		 covCombine() {};
  shared_obs(SEXP x) {
    using tmbutils::asArray;
    if(Rf_isNull(getListElement(x,"hasSharedObs")) ||
       (int)*REAL(getListElement(x,"hasSharedObs")) == 0){
      hasSharedObs = 0;
    }else{
      hasSharedObs = (int)*REAL(getListElement(x,"hasSharedObs"));
      fleetTypes = asVector<int>(getListElement(x,"fleetTypes"));
      maxAgePlusGroup = asVector<int>(getListElement(x,"maxAgePlusGroup"));
      aux = asArray<int>(getListElement(x,"aux"));   
      logobs = asVector<Type>(getListElement(x,"logobs"));
      keyFleetStock = asMatrix<Type>(getListElement(x,"keyFleetStock"));   
      noYears = (int)*REAL(getListElement(x,"noYears"));
      noFleets = (int)*REAL(getListElement(x,"noFleets"));
      idx1 = asArray<int>(getListElement(x,"idx1"));
      idx2 = asArray<int>(getListElement(x,"idx2"));
      idxCor = asArray<int>(getListElement(x,"idxCor"));
      vector<int> covCombineTmp = asVector<int>(getListElement(x,"covCombine"));
      covCombine = vector<CovCombineType>(covCombineTmp.size());
      for(int i = 0; i < covCombine.size(); ++i)
	covCombine(i) = static_cast<CovCombineType>(covCombineTmp(i));
    }
  };  
};
