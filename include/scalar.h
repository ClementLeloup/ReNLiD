// scalar.h

#ifndef _SCALAR_H_
#define _SCALAR_H_

//#include <cmath>
#include <complex>
#include <vector>
//#include <iomanip>
// using namespace std::complex_literals;

#include "precisions.h"
#include "field.h"

// Higgs field with fundamental representation
class FundHiggs : public Field<std::complex<FloatType>>, public OperationHelperField<FundHiggs,std::complex<FloatType>>
{
public:
  // Constructors
  FundHiggs(std::vector<std::complex<FloatType>> phi);
  FundHiggs(){}
  
  FundHiggs conj(); // Get conjugate of the field
  FloatType V(); // Get field potential
};


// Higgs field with adjoint representation
class AdHiggs : public Field<FloatType>, public OperationHelperField<AdHiggs,FloatType>
{
public:
  //Constructors
  AdHiggs(std::vector<FloatType> phi);
  AdHiggs(){}

  FloatType V(); // Get field potential
};


// Simple real field
class RealScalar : public Field<FloatType>, public OperationHelperField<RealScalar,FloatType>
{
public:
  //Constructors
  RealScalar(std::vector<FloatType> phi);
  RealScalar(){}

  FloatType V(); // Get field potential
};


#endif
