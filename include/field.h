// field.h 

#ifndef _FIELD_H_
#define _FIELD_H_

#include <stdio.h>
#include <iostream>
#include <typeinfo>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>

#include "precisions.h"


// Class to propagate basic operations to all field child classes using curiously recurring template pattern (CRTP)
template<class S, class T>
struct OperationHelperField
{
  S operator*(const std::complex<FloatType>& s){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> sfield;//(n);
    for(int i=0; i<n; i++){
      sfield.push_back(s*reference[i]);
    }
    return S(sfield);
  }
  S operator*(const FloatType& s){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> sfield;//(n);
    for(int i=0; i<n; i++){
      sfield.push_back(s*reference[i]);
    }
    return S(sfield);
  }

  S operator/(const std::complex<FloatType>& s){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> sfield;//(n);
    for(int i=0; i<n; i++){
      sfield.push_back(reference[i]/s);
    }
    return S(sfield);
  }
  S operator/(const FloatType& s){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> sfield;//(n);
    for(int i=0; i<n; i++){
      sfield.push_back(reference[i]/s);
    }
    return S(sfield);
  }
  
  S operator+(const std::complex<FloatType>& s){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> sfield;//(n);
    for(int i=0; i<n; i++){
      sfield.push_back(reference[i]+s);
    }
    return S(sfield);
  }
  S operator+(const FloatType& s){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> sfield;//(n);
    for(int i=0; i<n; i++){
      sfield.push_back(reference[i]+s);
    }
    return S(sfield);
  }
  
  S operator-(const std::complex<FloatType>& s){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> sfield;//(n);
    for(int i=0; i<n; i++){
      sfield.push_back(reference[i]-s);
    }
    return S(sfield);
  }
  S operator-(const FloatType& s){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> sfield;//(n);
    for(int i=0; i<n; i++){
      sfield.push_back(reference[i]-s);
    }
    return S(sfield);
  }

  T operator*(const S& phi2){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    T scalarproduct = 0;
    for(int i=0; i<n; i++){
      scalarproduct += reference[i] * phi2[i];
    }
    return scalarproduct;
  }

  S operator+(const S& phi2){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> sumfield;//(n);
    for(int i=0; i<n; i++){
      sumfield.push_back(reference[i] + phi2[i]);
    }
    return S(sumfield);
  }
  S operator-(const S& phi2){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> difffield;//(n);
    for(int i=0; i<n; i++){
      difffield.push_back(reference[i] - phi2[i]);
    }
    return S(difffield);
  }
  S operator-(){
    S& reference = static_cast<S&>(*this);
    int n = reference.getSize();
    std::vector<T> minusfield;//(n);
    for(int i=0; i<n; i++){
      minusfield.push_back(-reference[i]);
    }
    return S(minusfield);
  }

};


// Class of general fields (members of a vector space with scalar product)
template <class T>
class Field
{
public:  
  // Constructors
  Field(std::vector<T> phi) : fsize(phi.size()), field(phi) {}
  Field(uint n) : fsize(n), field(n) {}
  Field(){}

  // Get element i of the field
  T &operator[](int i);
  const T& operator[](int i) const;

  // Get dimension
  int getSize() const;

  // Compute norm of field from scalar product
  FloatType normSquared() const;
  FloatType norm() const;

  // Print useful information about the field
  void print();

protected:
  // Protected variables
  int fsize;
  std::vector<T> field;

};

#endif
