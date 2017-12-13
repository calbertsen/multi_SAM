

#ifndef _CMOE_PARAMETER_TYPES_
#define _CMOE_PARAMETER_TYPES_


template<class Type>
struct cmoe_vector{

  vector<vector<Type> > dat;
  int size_;
  int ncol_;

  cmoe_vector(SEXP x){
    SEXP dimsS = PROTECT(getAttrib(x,install("vdim")));
    int nv = Rf_length((dimsS));
    ncol_ = nv;
    dat = vector<vector<Type> >(ncol_);
    int indx = 0;
    for(int i = 0; i < nv; ++i){
      int n = (int)REAL(dimsS)[i];
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
  }

  cmoe_vector(vector<vector<Type> > x){
    dat = x;
    int nv = x.size();
    ncol_ = nv;
    int res = 0;
    for(int i = 0; i < nv; ++i)
      res += x(i).size();
    size_ = res;
  }

  int size(){
    return size_;
  }

  vector<Type> col(int i){
    return dat(i);
  }

  int cols(){
    return ncol_;
  }

  Type sum(){
    Type res = Type(0.0);
    for(int i = 0; i < size_; ++i)
      res += (*this)[i];
    return res;
  }

  Type& operator()(int i){
    return (*this)[i];
  }

  // This can probably be done in a better way!
  Type& operator[](int i){
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

  // Addition with scalars
  cmoe_vector<Type>& operator +=(const Type& rhs){
    for(int i = 0; i < (*this).size(); ++i)
      (*this)[i] += rhs;
    return *this;
  }
  friend cmoe_vector<Type> operator+(cmoe_vector<Type> lhs, const Type& rhs){
    lhs += rhs;
    return lhs;
  }
  // Subtraction with scalars
  cmoe_vector<Type>& operator -=(const Type& rhs){
    for(int i = 0; i < (*this).size(); ++i)
      (*this)[i] -= rhs;
    return *this;
  }
  friend cmoe_vector<Type> operator-(cmoe_vector<Type> lhs, const Type& rhs){
    lhs -= rhs;
    return lhs;
  }
  // Multiplication with scalars
  cmoe_vector<Type>& operator *=(const Type& rhs){
    for(int i = 0; i < (*this).size(); ++i)
      (*this)[i] *= rhs;
    return *this;
  }
  friend cmoe_vector<Type> operator*(cmoe_vector<Type> lhs, const Type& rhs){
    lhs *= rhs;
    return lhs;
  }
  // Division with scalars
  cmoe_vector<Type>& operator /=(const Type& rhs){
    for(int i = 0; i < (*this).size(); ++i)
      (*this)[i] /= rhs;
    return *this;
  }
  friend cmoe_vector<Type> operator/(cmoe_vector<Type> lhs, const Type& rhs){
    lhs /= rhs;
    return lhs;
  }


};

template<class Type>
cmoe_vector<Type> asCmoeVector(SEXP x){
  return cmoe_vector<Type>(x);
}

template<class Type>
SEXP asSEXP(const cmoe_vector<Type>& x){
  return asSEXP(x.dat);
}


template<class Type>
struct cmoe_matrix{

  vector<matrix<Type> > dat;
  int size_;
  int ncol_;

  cmoe_matrix(SEXP x){
    SEXP dimcS = PROTECT(getAttrib(x,install("cdim")));
    SEXP dimrS = PROTECT(getAttrib(x,install("rdim")));
    int nv = Rf_length(dimcS);
    ncol_ = nv;
    dat = vector<matrix<Type> >(ncol_);
    int indx = 0;
    for(int i = 0; i < nv; ++i){
      int nr = REAL(dimrS)[i];
      int nc = REAL(dimcS)[i];
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

  cmoe_matrix(vector<matrix<Type> > x){
    dat = x;
    int nv = x.size();
    ncol_ = nv;
    int res = 0;
    for(int i = 0; i < nv; ++i)
      res += x(i).size();
    size_ = res;
  }

  int size(){
    return size_;
  }

  matrix<Type> col(int i){
    return dat(i);
  }

  int cols(){
    return ncol_;
  }

  Type sum(){
    Type res = Type(0.0);
    for(int i = 0; i < size_; ++i)
      res += (*this)[i];
    return res;
  }

  Type& operator()(int i){
    return (*this)[i];
  }

  // This can probably be done in a better way!
  Type& operator[](int i){
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


 // Addition with scalars
  cmoe_matrix<Type>& operator +=(const Type& rhs){
    for(int i = 0; i < (*this).size(); ++i)
      (*this)[i] += rhs;
    return *this;
  }
  friend cmoe_matrix<Type> operator+(cmoe_matrix<Type> lhs, const Type& rhs){
    lhs += rhs;
    return lhs;
  }
  // Subtraction with scalars
  cmoe_matrix<Type>& operator -=(const Type& rhs){
    for(int i = 0; i < (*this).size(); ++i)
      (*this)[i] -= rhs;
    return *this;
  }
  friend cmoe_matrix<Type> operator-(cmoe_matrix<Type> lhs, const Type& rhs){
    lhs -= rhs;
    return lhs;
  }
  // Multiplication with scalars
  cmoe_matrix<Type>& operator *=(const Type& rhs){
    for(int i = 0; i < (*this).size(); ++i)
      (*this)[i] *= rhs;
    return *this;
  }
  friend cmoe_matrix<Type> operator*(cmoe_matrix<Type> lhs, const Type& rhs){
    lhs *= rhs;
    return lhs;
  }
  // Division with scalars
  cmoe_matrix<Type>& operator /=(const Type& rhs){
    for(int i = 0; i < (*this).size(); ++i)
      (*this)[i] /= rhs;
    return *this;
  }
  friend cmoe_matrix<Type> operator/(cmoe_matrix<Type> lhs, const Type& rhs){
    lhs /= rhs;
    return lhs;
  }


  
};

template<class Type>
cmoe_matrix<Type> asCmoeMatrix(SEXP x){
  return cmoe_matrix<Type>(x);
}

template<class Type>
SEXP asSEXP(const cmoe_matrix<Type>& x){
  return asSEXP(x.dat);
}


#define PARAMETER_CMOE_VECTOR(name) cmoe_vector<Type> name(objective_function::fillShape(asCmoeVector<Type>(objective_function::getShape(#name,&isArray)),#name));

#define PARAMETER_CMOE_MATRIX(name) cmoe_matrix<Type> name(objective_function::fillShape(asCmoeMatrix<Type>(objective_function::getShape(#name,&isArray)),#name));





template <class Type>
vector<vector<Type> > exp(cmoe_vector<Type> x){
  cmoe_vector<Type> tmp(x.dat);
  for(int i = 0; i < tmp.size(); ++i)
    tmp[i] = exp(tmp[i]);
  return tmp.dat;
}



#endif
