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
FundHiggs::FundHiggs(std::vector<std::complex<FloatType>> phi){
  if(phi.size()!=2){
    fprintf(stderr,"size is %ld\n", phi.size());
    throw std::invalid_argument("FundHiggs has the wrong size.");
  } else{
    field=phi;
    fsize=phi.size();
  }
}

// Get conjugate of the field
FundHiggs FundHiggs::conj(){
  std::vector<std::complex<FloatType>> conjfield;//(size);// = new std::complex<FloatType>[size];
  for(int i=0; i<fsize; i++){
    // conjfield[i] = std::conj(field[i]);
    conjfield.push_back(field[i]);
  }
  return FundHiggs(conjfield);
}



// Higgs field with adjoint representation
//Constructors
AdHiggs::AdHiggs(std::vector<FloatType> phi){
  if(phi.size()!=3){
    fprintf(stderr,"size is %ld\n", phi.size());
    throw std::invalid_argument("AddHiggs has the wrong size.");
  } else{
    field=phi;
    fsize=phi.size();
  }
}


// Simple real scalar field
//Constructors
RealScalar::RealScalar(std::vector<FloatType> phi){
  if(phi.size()!=1){
    fprintf(stderr,"size is %ld\n", phi.size());
    throw std::invalid_argument("RealScalar has the wrong size.");
  } else{
    field=phi;
    fsize=phi.size();
  }
}

//Potential
FloatType RealScalar::V(){
  RealScalar phi = (*this);
  return 0.25*(phi*phi)*(phi*phi);// (this*this)*(this*this);
  // return lamb*(phi*phi)*(phi*phi);
}

//Potential derivative
RealScalar RealScalar::dVdphi(){
  RealScalar phi = (*this);  
  return phi*(phi*phi); // this*(this*this);
  // return phi*4.0*lamb*(phi*phi);
}
