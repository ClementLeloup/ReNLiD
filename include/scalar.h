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
  // FundHiggs(std::complex<FloatType>* phi) : Field(phi, 2){}
  // FundHiggs(std::complex<FloatType>* phi, int n);
  FundHiggs(std::vector<std::complex<FloatType>> phi);
  FundHiggs(){}
  
  // Get conjugate of the field
  FundHiggs conj();
};




// Higgs field with adjoint representation
class AdHiggs : public Field<FloatType>, public OperationHelperField<AdHiggs,FloatType>
{
public:
  //Constructors
  // AdHiggs(FloatType* phi) : Field(phi, 3){}
  // AdHiggs(FloatType* phi, int n);
  AdHiggs(std::vector<FloatType> phi);
  AdHiggs(){}
};


#endif
