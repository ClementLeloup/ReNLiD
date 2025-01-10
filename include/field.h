// field.h 
#ifndef _FIELD_H_
#define _FIELD_H_

#include "precisions.h"

#include <stdio.h>
#include <iostream>
#include <typeinfo>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>


// Class to propagate basic operations to all field child classes
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

  // S& operator=(const S& phi2){
  //   S& reference = static_cast<S&>(*this);
  //   int n = reference.getSize();
  //   std::vector<T> copyfield;//(n);
  //   for(int i=0; i<n; i++){
  //     copyfield.push_back(phi2[i]);
  //   }
  //   return S(copyfield);
  // }
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


// Class of fields
template <class T>
class Field
{
public:  
  // Constructors
  // Field(T* phi, int n);
  // Field(T* phi, int n) : size(n), field(phi, phi+n){}
  Field(std::vector<T> phi) : fsize(phi.size()), field(phi) {}
  Field(uint n) : fsize(n), field(n) {}
  Field(){}

  // Get element i of the field
  T &operator[](int i);
  const T& operator[](int i) const;

  // Get dimension
  int getSize() const;

  // Compute norm of field
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
