// group.h

#ifndef _GROUP_H_
#define _GROUP_H_

#include<vector>

#include "precisions.h"
#include "field.h"
#include "lie.h"
#include "vector.h"
#include "scalar.h"

// U(1) group representation
class U1 : public Lie<u1>
{
public:
  // Constructors
  U1(u1 A) : Lie(A){}
  U1(FloatType theta) : Lie(u1(theta)){}
  U1(){}

  // Group multiplication
  std::complex<FloatType> action(const std::complex<FloatType>& a);
  U1 groupMultiply(const U1& U);
  U1 inverse();
};


// SU(2) group representation
class SU2 : public Lie<su2>
{
public:
  // Constructors
  SU2(su2 A) : Lie(A){}
  SU2(std::vector<FloatType> A) : Lie(su2(A)){}
  SU2(){}
  
  // Group action on fundamental representation
  FundHiggs action(const FundHiggs& h);

  // Group action on adjoint representation
  AdHiggs action(const AdHiggs& h);

  // Cross-product between SU2 algebra elements taken as SO3 representations
  SU2 crossProduct(const SU2& U);

  // Group inverse
  SU2 inverse();
  
  // Group multiplication
  SU2 groupMultiply(const SU2& U);

  // Print matrix form
  void matrixForm();

  void productMatrixForm(const SU2& U);
};


#endif
