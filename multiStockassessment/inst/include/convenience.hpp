template<class Type>
matrix<Type> cov2cor(matrix<Type> x){
  matrix<Type> Is(x.rows(),x.cols());
  Is.setZero();
  for(int i = 0; i < Is.rows(); ++i)
    Is(i,i) = Type(1.0)/sqrt(x(i,i));

  matrix<Type> r = (matrix<Type>)(Is * x) * Is;
  return r;
}



/* Norm of a vector */
template<class Type>
Type norm(vector<Type> x){
  Type res = 0.0;
  for(int i = 0; i < x.size(); ++i)
    res += x(i)*x(i);
  return sqrt(res);
}

/* Inner product of two vectors */
template<class Type>
Type vprod(vector<Type> x,vector<Type> y){
  Type res = 0.0;
  for(int i = 0; i < x.size(); ++i)
    res += x(i)*y(i);
  return res;
}
