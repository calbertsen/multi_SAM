MSAM_DEPENDS(convenience)

/* Structure to define covariance constraints */
HEADER(
template<class Type>
struct cov_constraints {
  vector<vector<Type> > a;
  cov_constraints(SEXP x);
};
       )

SOURCE(
       template<class Type>
       cov_constraints<Type>::cov_constraints(SEXP x){ // Constructor
	 a = vector<vector<Type> >(Rf_length(x));
	 for(int i = 0; i < a.size(); ++i){
	   vector<Type> tmp = asVector<Type>(VECTOR_ELT(x,i));
	   //a(i) = vector<Type>(tmp.size());
	   a(i) = tmp;
	 }
       }
       );
       /* Structure to define covariance design matrices */

MSAM_SPECIALIZATION(struct cov_constraints<double>);
MSAM_SPECIALIZATION(struct cov_constraints<TMBad::ad_aug>);


template<class Type>
matrix<Type> constructL(int nages, int nAreas, vector<Type> RE, matrix<Type> X, cov_constraints<Type> cons)SOURCE({

    int nAge = nages * nAreas;
    matrix<Type> llt(nAge,nAge);
    llt.setIdentity();
    int i,j,k=0;
    vector<Type> REX = X * RE;
    if(RE.size() > 0)
      for(i=0;i<llt.cols();i++){
	for(j=0;j<llt.rows();j++){
	  if(i<j){llt(j,i)=REX(k++);} // LAge becomes lower triangular
	}
      }

    for(int i = 0; i < llt.rows(); ++i)
      llt.row(i) /= norm((vector<Type>)llt.row(i));

    // Modify to get constraints
    matrix<Type> L = llt;
    for(int i = 0; i < cons.a.size(); ++i){
      matrix<Type> Y(cons.a(i).size()-1,nAge);
      for(int j = 1; j < cons.a(i).size(); ++j)
	Y.row(j-1) = L.row(CppAD::Integer(cons.a(i)(j)));
      
      for(int j = 0; j < Y.rows(); ++j){
	Y.row(j) /= norm((vector<Type>)Y.row(j));
	for(int k = j+1; k < Y.rows(); ++k){
	  vector<Type> tmp0 = (vector<Type>)Y.row(k) - vprod((vector<Type>)Y.row(k),((vector<Type>)Y.row(j))) * (vector<Type>)Y.row(j);
	  Y.row(k) = tmp0;
	}
      }
      vector<Type> tmp = L.row(CppAD::Integer(cons.a(i)(0)));
      for(int j = 0; j < Y.rows(); ++j){
	vector<Type> t0 = vprod(tmp,(vector<Type>)Y.row(j))*(vector<Type>)Y.row(j);
	tmp -= t0;
      }
      L.row(CppAD::Integer(cons.a(i)(0))) = tmp / norm(tmp);
    }

    return L;
  })

MSAM_SPECIALIZATION(matrix<double> constructL(int, int, vector<double>, matrix<double>, cov_constraints<double>));
MSAM_SPECIALIZATION(matrix<TMBad::ad_aug> constructL(int, int, vector<TMBad::ad_aug>, matrix<TMBad::ad_aug>, cov_constraints<TMBad::ad_aug>));
