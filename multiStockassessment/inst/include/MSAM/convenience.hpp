MSAM_DEPENDS(param_types)

template<class Type>
matrix<Type> cov2cor(matrix<Type> x)SOURCE({
  matrix<Type> Is(x.rows(),x.cols());
  Is.setZero();
  for(int i = 0; i < Is.rows(); ++i)
    Is(i,i) = Type(1.0)/sqrt(x(i,i));

  matrix<Type> r = (matrix<Type>)(Is * x) * Is;
  return r;
  });

MSAM_SPECIALIZATION(matrix<double> cov2cor(matrix<double>));
MSAM_SPECIALIZATION(matrix<TMBad::ad_aug> cov2cor(matrix<TMBad::ad_aug>));



/* Norm of a vector */
template<class Type>
Type norm(vector<Type> x)SOURCE({
  Type res = 0.0;
  for(int i = 0; i < x.size(); ++i)
    res += x(i)*x(i);
  return sqrt(res);
  });

MSAM_SPECIALIZATION(double norm(vector<double>));
MSAM_SPECIALIZATION(TMBad::ad_aug norm(vector<TMBad::ad_aug>));


/* Inner product of two vectors */
template<class Type>
Type vprod(vector<Type> x,vector<Type> y)SOURCE({
  Type res = 0.0;
  for(int i = 0; i < x.size(); ++i)
    res += x(i)*y(i);
  return res;
  });

MSAM_SPECIALIZATION(double vprod(vector<double>, vector<double>));
MSAM_SPECIALIZATION(TMBad::ad_aug vprod(vector<TMBad::ad_aug>, vector<TMBad::ad_aug>));



template<class Type>
array<Type> getArray(cmoe_matrix<Type> x, int s)SOURCE({
  array<Type> Xa(x.col(s).rows(),x.col(s).cols());
  Xa = x.col(s);
  return Xa;
  })

MSAM_SPECIALIZATION(array<double> getArray(cmoe_matrix<double>, int));
MSAM_SPECIALIZATION(array<TMBad::ad_aug> getArray(cmoe_matrix<TMBad::ad_aug>, int));

template<class Type>
array<Type> getArray(cmoe_3darray<Type> x, int s)SOURCE({
  // vector<int> dim = x.col(s).dim;
  // array<Type> Xa(dim[0], dim[1], dim[2]);
  // Xa = x.col(s);
  // return Xa;
  return x.col(s);
  })

MSAM_SPECIALIZATION(array<double> getArray(cmoe_3darray<double>, int));
MSAM_SPECIALIZATION(array<TMBad::ad_aug> getArray(cmoe_3darray<TMBad::ad_aug>, int));
