
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
  int hasSharedObs;
  vector<int> fleetTypes;
  vector<int> maxAgePlusGroup;
  array<int> aux;
  vector<Type> logobs;
  matrix<Type> keyFleetStock;

  shared_obs() : hasSharedObs(0),
		 fleetTypes(),
		 maxAgePlusGroup(),
		 aux(),
		 logobs(),
		 keyFleetStock() {};
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
    }
  };
};
