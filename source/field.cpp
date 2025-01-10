// field.cpp

#include "field.h"
#include "precisions.h"
#include <complex>
#include <vector>

// Get element i of the field
template<class T>
T& Field<T>::operator[](int i){
  return field[i];
}

template<class T>
const T& Field<T>::operator[](int i) const {
  return field[i];
}

template<class T>
int Field<T>::getSize() const {
  return fsize;
}

template<class T>
FloatType Field<T>::normSquared() const {
  FloatType normsquared = 0;
  for(int i=0; i<fsize; i++){
    normsquared += pow(std::abs(field[i]), 2);
  }
  return normsquared;
}

template<class T>
FloatType Field<T>::norm() const {
  return sqrt(normSquared());
}
  
template<class T>
void Field<T>::print(){
  std::cout << "dimension = " << fsize << "\t";
  std::cout << "field = ";
  for(int i=0; i<fsize; i++){
    std::cout << "   " << field[i];
  }
  std::cout << "\n";
}


// Explicit instantiations
template class Field<FloatType>;
template class Field<std::complex<FloatType>>;
// You will only be able to use Field with real or complex
