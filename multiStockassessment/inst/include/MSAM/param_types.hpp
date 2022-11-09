
HEADER(
       template<class Type>
       struct cmoe_vector{

	 vector<vector<Type> > dat;
	 int size_;
	 int ncol_;

	 cmoe_vector(SEXP x);
	 cmoe_vector(vector<vector<Type> > x);

	 int size();

	 vector<Type>& col(int i);

	 int cols();

	 Type sum();

	 Type& operator()(int i);

	 // This can probably be done in a better way!
	 Type& operator[](int i);

	 // Addition with scalars
	 cmoe_vector<Type>& operator +=(const Type& rhs);
	 cmoe_vector<Type> operator+(const Type& rhs);
	 // Subtraction with scalars
	 cmoe_vector<Type>& operator -=(const Type& rhs);
	 cmoe_vector<Type> operator-(const Type& rhs);
	 // Multiplication with scalars
	 cmoe_vector<Type>& operator *=(const Type& rhs);
	 cmoe_vector<Type> operator*(const Type& rhs);
	 // Division with scalars
	 cmoe_vector<Type>& operator /=(const Type& rhs);
	 cmoe_vector<Type> operator/(const Type& rhs);

       };
       )


SOURCE(
       template<class Type>
       cmoe_vector<Type>::cmoe_vector(SEXP x){
	 SEXP dimsS = PROTECT(getAttrib(x,install("vdim")));
	 int nv = Rf_length((dimsS));
	 ncol_ = nv;
	 dat = vector<vector<Type> >(ncol_);
	 int indx = 0;
	 for(int i = 0; i < nv; ++i){
	   int n = INTEGER(dimsS)[i];
	   //) = vector<Type>(n);
	   vector<Type> tmp(n);
	   for(int j = 0; j < n; ++j){
	     tmp(j) = REAL(x)[indx];
	     indx++;
	   }
	   dat(i) = tmp;
      
	 }
	 size_ = indx;
	 UNPROTECT(1);
       });

SOURCE(
       template<class Type>
       cmoe_vector<Type>::cmoe_vector(vector<vector<Type> > x){
	 dat = x;
	 int nv = x.size();
	 ncol_ = nv;
	 int res = 0;
	 for(int i = 0; i < nv; ++i)
	   res += x(i).size();
	 size_ = res;
       });

SOURCE(
       template<class Type>
       int cmoe_vector<Type>::size(){
	 return size_;
       });

SOURCE(
       template<class Type>
       vector<Type>& cmoe_vector<Type>::col(int i){
	 return dat(i);
       });

SOURCE(
       template<class Type>
       int cmoe_vector<Type>::cols(){
	 return ncol_;
       });

SOURCE(
       template<class Type>
       Type cmoe_vector<Type>::sum(){
	 Type res = Type(0.0);
	 for(int i = 0; i < size_; ++i)
	   res += (*this)[i];
	 return res;
       });

SOURCE(
       template<class Type>
       Type& cmoe_vector<Type>::operator()(int i){
	 return (*this)[i];
       }
       )

// This can probably be done in a better way!
SOURCE(
       template<class Type>
       Type& cmoe_vector<Type>::operator[](int i){
	 int dim;
	 for(int q = 0; q < ncol_; ++q){
	   dim = dat(q).size();
	   if(i < dim){
	     return(dat(q)(i));
	   }
	   i -= dim;
	 }
	 return dat(-1)(0);		// Not a pretty solution
       });

// Addition with scalars
SOURCE(
       template<class Type>
       cmoe_vector<Type>& cmoe_vector<Type>::operator +=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] += rhs;
	 return *this;
       });
SOURCE(
       template<class Type>
       cmoe_vector<Type> cmoe_vector<Type>::operator+(const Type& rhs){
	 cmoe_vector<Type> lhs(this->dat);
	 lhs += rhs;
	 return lhs;
       });
// Subtraction with scalars
SOURCE(
       template<class Type>
       cmoe_vector<Type>& cmoe_vector<Type>::operator -=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] -= rhs;
	 return *this;
       });

SOURCE(
       template<class Type>
       cmoe_vector<Type> cmoe_vector<Type>::operator-(const Type& rhs){
	 cmoe_vector<Type> lhs(this->dat);
	 lhs -= rhs;
	 return lhs;
       });

// Multiplication with scalars
SOURCE(
       template<class Type>
       cmoe_vector<Type>& cmoe_vector<Type>::operator *=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] *= rhs;
	 return *this;
       }
       );

SOURCE(
       template<class Type>
       cmoe_vector<Type> cmoe_vector<Type>::operator*(const Type& rhs){
	 cmoe_vector<Type> lhs(this->dat);
	 lhs *= rhs;
	 return lhs;
       });

// Division with scalars
SOURCE(
       template<class Type>
       cmoe_vector<Type>& cmoe_vector<Type>::operator /=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] /= rhs;
	 return *this;
       });

SOURCE(
       template<class Type>
       cmoe_vector<Type> cmoe_vector<Type>::operator/(const Type& rhs){
	 cmoe_vector<Type> lhs(this->dat);
	 lhs /= rhs;
	 return lhs;
       });
  


MSAM_SPECIALIZATION(struct cmoe_vector<double>);
MSAM_SPECIALIZATION(struct cmoe_vector<TMBad::ad_aug>);

template<class Type>
cmoe_vector<Type> asCmoeVector(SEXP x)SOURCE({
    return cmoe_vector<Type>(x);
  });

MSAM_SPECIALIZATION(cmoe_vector<double> asCmoeVector(SEXP));
MSAM_SPECIALIZATION(cmoe_vector<TMBad::ad_aug> asCmoeVector(SEXP));


template<class Type>
SEXP asSEXP(const cmoe_vector<Type>& x)SOURCE({
    return asSEXP(x.dat);
  });

MSAM_SPECIALIZATION(SEXP asSEXP(const cmoe_vector<double>&));
MSAM_SPECIALIZATION(SEXP asSEXP(const cmoe_vector<TMBad::ad_aug>&));


HEADER(
       template<class Type>
       struct cmoe_matrix{

	 vector<matrix<Type> > dat;
	 int size_;
	 int ncol_;

	 cmoe_matrix(SEXP x);
	 cmoe_matrix(vector<matrix<Type> > x);
	 cmoe_matrix(const cmoe_matrix<Type>& x);
  
	 int size() const;

	 matrix<Type>& col(int i);

	 int cols() const;

	 Type sum();

	 Type& operator()(int i);

	 // This can probably be done in a better way!
	 Type& operator[](int i);

	 // Addition with scalars
	 cmoe_matrix<Type>& operator +=(const Type& rhs);
	 cmoe_matrix<Type> operator+(const Type& rhs);
	 // Subtraction with scalars
	 cmoe_matrix<Type>& operator -=(const Type& rhs);
	 cmoe_matrix<Type> operator-(const Type& rhs);
	 // Multiplication with scalars
	 cmoe_matrix<Type>& operator *=(const Type& rhs);
	 cmoe_matrix<Type> operator*(const Type& rhs);
	 // Division with scalars
	 cmoe_matrix<Type>& operator /=(const Type& rhs);
	 cmoe_matrix<Type> operator/(const Type& rhs);
       };
       )

SOURCE(
       template<class Type>
       cmoe_matrix<Type>::cmoe_matrix(SEXP x){
	 SEXP dimcS = PROTECT(getAttrib(x,install("cdim")));
	 SEXP dimrS = PROTECT(getAttrib(x,install("rdim")));
	 int nv = Rf_length(dimcS);
	 ncol_ = nv;
	 dat = vector<matrix<Type> >(ncol_);
	 int indx = 0;
	 for(int i = 0; i < nv; ++i){
	   int nr = INTEGER(dimrS)[i];
	   int nc = INTEGER(dimcS)[i];
	   //) = vector<Type>(n);
	   matrix<Type> tmp(nr,nc);
	   for(int k = 0; k < nc; ++k)
	     for(int j = 0; j < nr; ++j){
	       tmp(j,k) = REAL(x)[indx];
	       indx++;
	     }
	   dat(i) = tmp;
      
	 }
	 size_ = indx;
	 UNPROTECT(2);
       }
       );

SOURCE(
       template<class Type>
       cmoe_matrix<Type>::cmoe_matrix(vector<matrix<Type> > x){    
	 dat = x;
	 int nv = x.size();
	 ncol_ = nv;
	 int res = 0;
	 for(int i = 0; i < nv; ++i)
	   res += x(i).size();
	 size_ = res;
       }
       );

SOURCE(
       template<class Type>
       cmoe_matrix<Type>::cmoe_matrix(const cmoe_matrix<Type>& x) : dat(x.dat),
       size_(x.size()),
       ncol_(x.cols()) {}
       );

SOURCE(
       template<class Type>
       int cmoe_matrix<Type>::size() const{
	 return size_;
       }
       );

SOURCE(
       template<class Type>
       matrix<Type>& cmoe_matrix<Type>::col(int i){
	 return dat(i);
       }
       );

SOURCE(
       template<class Type>
       int cmoe_matrix<Type>::cols() const{
	 return ncol_;
       }
       );

SOURCE(
       template<class Type>
       Type cmoe_matrix<Type>::sum(){
	 Type res = Type(0.0);
	 for(int i = 0; i < size_; ++i)
	   res += (*this)[i];
	 return res;
       }
       );

SOURCE(
       template<class Type>
       Type& cmoe_matrix<Type>::operator()(int i){
	 return (*this)[i];
       }
       );

SOURCE(
       template<class Type>
       // This can probably be done in a better way!
       Type& cmoe_matrix<Type>::operator[](int i){
	 int dim;
	 for(int q = 0; q < ncol_; ++q){
	   dim = dat(q).size();
	   if(i < dim){
	     return(dat(q)(i));
	   }
	   i -= dim;
	 }
	 return dat(-1)(0);		// Not a pretty solution
       }
       );


// Addition with scalars
SOURCE(
       template<class Type>
       cmoe_matrix<Type>& cmoe_matrix<Type>::operator +=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] += rhs;
	 return *this;
       }
       );


SOURCE(
       template<class Type>
       cmoe_matrix<Type> cmoe_matrix<Type>::operator+(const Type& rhs){
	 cmoe_matrix<Type> lhs(this->dat);
	 lhs += rhs;
	 return lhs;
       }
       );


// Subtraction with scalars
SOURCE(
       template<class Type>
       cmoe_matrix<Type>& cmoe_matrix<Type>::operator -=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] -= rhs;
	 return *this;
       }
       );


SOURCE(
       template<class Type>
       cmoe_matrix<Type> cmoe_matrix<Type>::operator-(const Type& rhs){
	 cmoe_matrix<Type> lhs(this->dat);
	 lhs -= rhs;
	 return lhs;
       }
       );

// Multiplication with scalars
SOURCE(
       template<class Type>
       cmoe_matrix<Type>& cmoe_matrix<Type>::operator *=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] *= rhs;
	 return *this;
       }
       );

SOURCE(
       template<class Type>
       cmoe_matrix<Type> cmoe_matrix<Type>::operator*(const Type& rhs){
	 cmoe_matrix<Type> lhs(this->dat);
	 lhs *= rhs;
	 return lhs;
       }
       );

// Division with scalars
SOURCE(
       template<class Type>
       cmoe_matrix<Type>& cmoe_matrix<Type>::operator /=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] /= rhs;
	 return *this;
       }
       );


SOURCE(
       template<class Type>
       cmoe_matrix<Type> cmoe_matrix<Type>::operator/(const Type& rhs){
	 cmoe_matrix<Type> lhs(this->dat);
	 lhs /= rhs;
	 return lhs;
       }  
       );


MSAM_SPECIALIZATION(struct cmoe_matrix<double>);
MSAM_SPECIALIZATION(struct cmoe_matrix<TMBad::ad_aug>);



template<class Type>
cmoe_matrix<Type> asCmoeMatrix(SEXP x)SOURCE({
    return cmoe_matrix<Type>(x);
  });

MSAM_SPECIALIZATION(cmoe_matrix<double> asCmoeMatrix(SEXP));
MSAM_SPECIALIZATION(cmoe_matrix<TMBad::ad_aug> asCmoeMatrix(SEXP));


template<class Type>
SEXP asSEXP(const cmoe_matrix<Type>& x)SOURCE({
    return asSEXP(x.dat);
  });


MSAM_SPECIALIZATION(SEXP asSEXP(const cmoe_matrix<double>&));
MSAM_SPECIALIZATION(SEXP asSEXP(const cmoe_matrix<TMBad::ad_aug>&));


// 3D Array
HEADER(
template<class Type>
struct cmoe_3darray{

  vector<array<Type> > dat;
  int size_;
  int ncol_;

  cmoe_3darray(SEXP x);
  cmoe_3darray(vector<array<Type> > x);

  template<class T>
  cmoe_3darray(cmoe_3darray<T> x);

  int size();
  array<Type>& col(int i);
  int cols();
  Type sum();

  Type& operator()(int i);

  // This can probably be done in a better way!
  Type& operator[](int i);


  // Addition with scalars
  cmoe_3darray<Type>& operator +=(const Type& rhs);
  cmoe_3darray<Type> operator+(const Type& rhs);
  // Subtraction with scalars
  cmoe_3darray<Type>& operator -=(const Type& rhs);
  cmoe_3darray<Type> operator-(const Type& rhs);
  // Multiplication with scalars
  cmoe_3darray<Type>& operator *=(const Type& rhs);
  cmoe_3darray<Type> operator*(const Type& rhs);
  // Division with scalars
  cmoe_3darray<Type>& operator /=(const Type& rhs);
  cmoe_3darray<Type> operator/(const Type& rhs);
  
};
       )


SOURCE(
       template<class Type>
       cmoe_3darray<Type>::cmoe_3darray(SEXP x){
	 SEXP dimcS = PROTECT(getAttrib(x,install("cdim")));
	 SEXP dimrS = PROTECT(getAttrib(x,install("rdim")));
	 SEXP dimaS = PROTECT(getAttrib(x,install("adim")));
	 int nv = Rf_length(dimcS);
	 ncol_ = nv;
	 dat = vector<array<Type> >(ncol_);
	 int indx = 0;
	 for(int i = 0; i < nv; ++i){
	   int nr = INTEGER(dimrS)[i];
	   int nc = INTEGER(dimcS)[i];
	   int na = INTEGER(dimaS)[i];
	   //) = vector<Type>(n);
	   array<Type> tmp(nr,nc,na);
	   tmp.setZero();
	   for(int a = 0; a < na; ++a)
	     for(int k = 0; k < nc; ++k)
	       for(int j = 0; j < nr; ++j){
		 tmp(j,k,a) = REAL(x)[indx];
		 indx++;
	       }
	   dat(i) = tmp;
      
	 }
	 size_ = indx;
	 UNPROTECT(3);
       }
       )


SOURCE(
       template<class Type>
       cmoe_3darray<Type>::cmoe_3darray(vector<array<Type> > x){
	 dat = x;
	 int nv = x.size();
	 ncol_ = nv;
	 int res = 0;
	 for(int i = 0; i < nv; ++i)
	   res += x(i).size();
	 size_ = res;
       }
       )


SOURCE(
       template<class Type>
       int cmoe_3darray<Type>::size(){
	 return size_;
       }
       )

SOURCE(
       template<class Type>
       array<Type>& cmoe_3darray<Type>::col(int i){
	 return dat(i);
       }
       )


SOURCE(
       template<class Type>
       int cmoe_3darray<Type>::cols(){
	 return ncol_;
       }
       )

SOURCE(
       template<class Type>
       Type cmoe_3darray<Type>::sum(){
	 Type res = Type(0.0);
	 for(int i = 0; i < size_; ++i)
	   res += (*this)[i];
	 return res;
       }
       )


SOURCE(
       template<class Type>
       Type& cmoe_3darray<Type>::operator()(int i){
	 return (*this)[i];
       }
       )


// This can probably be done in a better way!
SOURCE(
       template<class Type>
       Type& cmoe_3darray<Type>::operator[](int i){
	 int dim;
	 for(int q = 0; q < ncol_; ++q){
	   dim = dat(q).size();
	   if(i < dim){
	     return(dat(q)(i));
	   }
	   i -= dim;
	 }
	 return dat(-1)(0);		// Not a pretty solution
       }
       )

// Addition with scalars
SOURCE(
       template<class Type>
       cmoe_3darray<Type>& cmoe_3darray<Type>::operator +=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] += rhs;
	 return *this;
       }
       )


SOURCE(
       template<class Type>
       cmoe_3darray<Type> cmoe_3darray<Type>::operator+(const Type& rhs){
	 cmoe_3darray<Type> lhs(this->dat);
	 lhs += rhs;
	 return lhs;
       }
       )

// Subtraction with scalars
SOURCE(
       template<class Type>
       cmoe_3darray<Type>& cmoe_3darray<Type>::operator -=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] -= rhs;
	 return *this;
       }
       );

SOURCE(
       template<class Type>
       cmoe_3darray<Type> cmoe_3darray<Type>::operator-(const Type& rhs){
	 cmoe_3darray<Type> lhs(this->dat);
	 lhs -= rhs;
	 return lhs;
       }
       );

// Multiplication with scalars
SOURCE(
       template<class Type>
       cmoe_3darray<Type>& cmoe_3darray<Type>::operator *=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] *= rhs;
	 return *this;
       }
       );

SOURCE(
       template<class Type>
       cmoe_3darray<Type> cmoe_3darray<Type>::operator*(const Type& rhs){
	 cmoe_3darray<Type> lhs(this->dat);
	 lhs *= rhs;
	 return lhs;
       }
       );

// Division with scalars
SOURCE(
       template<class Type>
       cmoe_3darray<Type>& cmoe_3darray<Type>::operator /=(const Type& rhs){
	 for(int i = 0; i < (*this).size(); ++i)
	   (*this)[i] /= rhs;
	 return *this;
       }
       );

SOURCE(
       template<class Type>
       cmoe_3darray<Type> cmoe_3darray<Type>::operator/(const Type& rhs){
	 cmoe_3darray<Type> lhs(this->dat);
	 lhs /= rhs;
	 return lhs;
       }
       );  



MSAM_SPECIALIZATION(struct cmoe_3darray<double>);
MSAM_SPECIALIZATION(struct cmoe_3darray<TMBad::ad_aug>);


template<class Type>
cmoe_3darray<Type> asCmoe3darray(SEXP x)SOURCE({
  return cmoe_3darray<Type>(x);
  })


MSAM_SPECIALIZATION(cmoe_3darray<double> asCmoe3darray(SEXP));
MSAM_SPECIALIZATION(cmoe_3darray<TMBad::ad_aug> asCmoe3darray(SEXP));


// SEXP asSEXP(const cmoe_3darray<double>& x)SOURCE({
//     return asSEXP(x.dat);
//   })

//   SEXP asSEXP(const cmoe_3darray<TMBad::ad_aug>& x)SOURCE({
//     return asSEXP(x.dat);
//   })


// MSAM_SPECIALIZATION(SEXP asSEXP(const cmoe_3darray<double>&));
// MSAM_SPECIALIZATION(SEXP asSEXP(const cmoe_3darray<TMBad::ad_aug>&));




#define PARAMETER_CMOE_VECTOR(name) cmoe_vector<Type> name(objective_function::fillShape(asCmoeVector<Type>(objective_function::getShape(#name,&isArray)),#name));

#define PARAMETER_CMOE_MATRIX(name) cmoe_matrix<Type> name(objective_function::fillShape(asCmoeMatrix<Type>(objective_function::getShape(#name,&isArray)),#name));

#define PARAMETER_CMOE_3DARRAY(name) cmoe_3darray<Type> name(objective_function::fillShape(asCmoe3darray<Type>(objective_function::getShape(#name,&isArray)),#name));





template <class Type>
vector<vector<Type> > exp(cmoe_vector<Type> x)SOURCE({
  cmoe_vector<Type> tmp(x.dat);
  for(int i = 0; i < tmp.size(); ++i)
    tmp[i] = exp(tmp[i]);
  return tmp.dat;
  })


MSAM_SPECIALIZATION(vector<vector<double> > exp(const cmoe_vector<double>));
MSAM_SPECIALIZATION(vector<vector<TMBad::ad_aug> > exp(const cmoe_vector<TMBad::ad_aug>));
