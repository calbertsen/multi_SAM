SAM_DEPENDS(define)
SAM_DEPENDS(refpointset)
SAM_DEPENDS(forecast)

// Define data structure

HEADER(
template<class Type>
struct sam_data {

  vector<dataSet<Type> > dataSets;
  vector<confSet> confSets;
  vector<forecastSet<Type> > forecastSets;
  vector<referencepointList<Type> > referencepointLists;
  
  sam_data(SEXP x); // End constructor

};
       )

SOURCE(
       template<class Type>
       sam_data<Type>::sam_data(SEXP x) : dataSets(Rf_length(x)), confSets(Rf_length(x)),
       forecastSets(Rf_length(x)), referencepointLists(Rf_length(x)){ // Constructor
	 int n = Rf_length(x);

	 // Loop through x to get the data for each area
	 for(int i = 0; i < n; ++i){
	   SEXP y = PROTECT(VECTOR_ELT(x,i));
	   dataSets(i) = dataSet<Type>(y);
	   confSets(i) = confSet(y);      
	   forecastSets(i) = forecastSet<Type>(getListElement(y,"forecast"));
	   referencepointLists(i) = referencepointList<Type>(getListElement(y,"referencepoints"));
	   UNPROTECT(1);
	 }
       };
       )

MSM_SPECIALIZATION(struct sam_data<double>);
MSM_SPECIALIZATION(struct sam_data<TMBad::ad_aug>);


HEADER(
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
  array<Type> auxData;
  vector<Type> logobs;
  vector<Type> weight;
  matrix<Type> keyFleetStock;
  int noYears;
  int noFleets;
  array<int> idx1;
  array<int> idx2;
  array<int> idxCor;
  listMatrixFromR<Type> corList;

  

  vector<CovCombineType> covCombine;

  shared_obs();
  shared_obs(SEXP x);  
};
       )

SOURCE(
       template<class Type>
       shared_obs<Type>::shared_obs() : hasSharedObs(0),
       fleetTypes(),
       maxAgePlusGroup(),
       aux(),
       auxData(),
       logobs(),
       weight(),
       keyFleetStock(),
       noYears(0),
       noFleets(0),
       idx1(),
       idx2(),
       idxCor(),
       corList(),
       covCombine() {};
       )


SOURCE(
       template<class Type>
       shared_obs<Type>::shared_obs(SEXP x) {
	 using tmbutils::asArray;
	 if(Rf_isNull(getListElement(x,"hasSharedObs")) ||
	    (int)*REAL(getListElement(x,"hasSharedObs",&Rf_isNumeric)) == 0){
	   hasSharedObs = 0;
	 }else{
	   hasSharedObs = (int)*REAL(getListElement(x,"hasSharedObs",&Rf_isNumeric));
	   fleetTypes = asVector<int>(getListElement(x,"fleetTypes",&Rf_isNumeric));
	   maxAgePlusGroup = asVector<int>(getListElement(x,"maxAgePlusGroup",&Rf_isNumeric));
	   aux = asArray<int>(getListElement(x,"aux",&Rf_isArray));
	   auxData = asArray<Type>(getListElement(x,"auxData",&Rf_isArray));   
	   logobs = asVector<Type>(getListElement(x,"logobs",&Rf_isNumeric));
	   weight = asVector<Type>(getListElement(x,"weight",&Rf_isNumeric));
	   keyFleetStock = asMatrix<Type>(getListElement(x,"keyFleetStock",&Rf_isMatrix));   
	   noYears = (int)*REAL(getListElement(x,"noYears",&isNumericScalar));
	   noFleets = (int)*REAL(getListElement(x,"noFleets",&isNumericScalar));
	   idx1 = asArray<int>(getListElement(x,"idx1",&Rf_isArray));
	   idx2 = asArray<int>(getListElement(x,"idx2",&Rf_isArray));
	   idxCor = asArray<int>(getListElement(x,"idxCor", &Rf_isArray));
	   vector<int> covCombineTmp = asVector<int>(getListElement(x,"covCombine", &Rf_isNumeric));
	   corList = listMatrixFromR<Type>(getListElement(x,"corList"));
	   covCombine = vector<CovCombineType>(covCombineTmp.size());
	   for(int i = 0; i < covCombine.size(); ++i)
	     covCombine(i) = static_cast<CovCombineType>(covCombineTmp(i));
	 }
       };  
       )


MSM_SPECIALIZATION(struct shared_obs<double>);
MSM_SPECIALIZATION(struct shared_obs<TMBad::ad_aug>);
