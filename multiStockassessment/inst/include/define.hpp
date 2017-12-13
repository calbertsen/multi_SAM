
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



