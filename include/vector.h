// vector.h

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include<vector>

#include "precisions.h"
#include "field.h"


// U(1) algebra representation
class u1: public Field<FloatType>, public OperationHelperField<u1,FloatType>
{
public:
  u1(FloatType theta);
  u1(std::vector<FloatType> theta);
  u1(){}
};


// SU(2) algebra representation
class su2: public Field<FloatType>, public OperationHelperField<su2,FloatType>
{
public:
  su2(std::vector<FloatType> A);
  su2(){}
};

#endif
