// scalar.cpp

#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
// using namespace std::complex_literals;

#include "precisions.h"
#include "scalar.h"


// Higgs field with fundamental representation
// Constructors
// FundHiggs::FundHiggs(std::complex<FloatType>* phi, int n){
//   if(n!=2){
//     throw std::invalid_argument("FundHiggs has the wrong size.");
//   } else{
//     field=std::vector<std::complex<FloatType>> (phi, phi+n);
//     size=n;
//   }
// }

FundHiggs::FundHiggs(std::vector<std::complex<FloatType>> phi){

  
  if(phi.size()!=2){
    fprintf(stderr,"size is %ld\n", phi.size());
    throw std::invalid_argument("FundHiggs has the wrong size.");
  } else{
    field=phi;
    size=phi.size();
  }
}

// Get conjugate of the field
FundHiggs FundHiggs::conj(){
  std::vector<std::complex<FloatType>> conjfield;//(size);// = new std::complex<FloatType>[size];
  for(int i=0; i<size; i++){
    // conjfield[i] = std::conj(field[i]);
    conjfield.push_back(field[i]);
  }
  return FundHiggs(conjfield);
}



// Higgs field with adjoint representation
//Constructors
// AdHiggs::AdHiggs(FloatType* phi, int n){
//   if(n!=3){
//     throw std::invalid_argument("AdHiggs has the wrong size.");
//   } else{
//     field=phi;
//     size=n;
//   }
// }

AdHiggs::AdHiggs(std::vector<FloatType> phi){
  if(phi.size()!=3){
    fprintf(stderr,"size is %ld\n", phi.size());
    throw std::invalid_argument("AddHiggs has the wrong size.");
  } else{
    field=phi;
    size=phi.size();
  }
}
