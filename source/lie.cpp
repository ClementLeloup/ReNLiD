// lie.cpp

#include <assert.h>
using namespace std;

#include "field.h"
#include "precisions.h"
#include "lie.h"
#include "scalar.h"
#include "vector.h"


template<class T>
T Lie<T>::toAlgebra() const {
  return alg;
}

template<class T>
FloatType Lie<T>::algebraNorm() const {
  return alg.norm();
}

// To be cleaned
template<class T>
Lie<T>& Lie<T>::inverse(){

  assert(1==0); // group multiplication should be implemented at the child class level
  
  return (*this);
}

// To be cleaned
template<class T>
Lie<T>& Lie<T>::groupMultiply(const Lie<T>& U){

  assert(1==0); // group multiplication should be implemented at the child class level
  
  return (*this);
}

template<class T>
void Lie<T>::print(){
  std::cout << "Algebra element : \n";
  alg.print();
}


// Explicit instantiations
template class Lie<u1>;
template class Lie<su2>;
// You will only be able to use Field with real or complex
