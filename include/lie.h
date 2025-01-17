// lie.h

#ifndef _LIE_H_
#define _LIE_H_

#include <vector>

#include "field.h"
#include "precisions.h"

// Class of Lie group representation
// It is used to define the edges
template<class T>
class Lie
{
public:
  // Constructors
  Lie(T A){alg = A;};
  Lie(std::vector<FloatType> A) : alg(A){}
  Lie(){}

  // Get Lie algebra representation
  T toAlgebra() const;

  // Get norm
  FloatType algebraNorm() const;

  // Get conjugate
  Lie& inverse();

  // Group multiplication
  Lie& groupMultiply(const Lie& U); // this needs to be implemented in child classes
  
  // Print useful information about the field
  void print();

protected:
  // Protected variables
  T alg;

};


#endif
