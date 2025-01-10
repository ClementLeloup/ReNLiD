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

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>



// Higgs field with fundamental representation
class FundHiggs : public Field<std::complex<FloatType>>, public OperationHelperField<FundHiggs,std::complex<FloatType>>
{
public:
  // Constructors
  FundHiggs(std::vector<std::complex<FloatType>> phi);
  FundHiggs() : Field<std::complex<FloatType>>(2) {}
  
  FundHiggs conj(); // Get conjugate of the field
  FloatType V(); // Get field potential
  FloatType dVdphi(); // Get potential derivative wrt field
};


// Higgs field with adjoint representation
class AdHiggs : public Field<FloatType>, public OperationHelperField<AdHiggs,FloatType>
{
public:
  //Constructors
  AdHiggs(std::vector<FloatType> phi);
  AdHiggs() : Field<FloatType>(3) {}

  FloatType V(); // Get field potential
  FloatType dVdphi(); // Get potential derivative wrt field
};


// Simple real field
class RealScalar : public Field<FloatType>, public OperationHelperField<RealScalar,FloatType>
{

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & field;
    ar & fsize;
  }
  
public:
  //Constructors
  RealScalar(std::vector<FloatType> phi);
  RealScalar() : Field<FloatType>(1) {}

  FloatType V(); // Get field potential
  RealScalar dVdphi(); // Get potential derivative wrt field
};


#endif
