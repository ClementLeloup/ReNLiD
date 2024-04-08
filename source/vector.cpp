// vector.cpp

#include <vector>

#include "precisions.h"
#include "scalar.h"
#include "vector.h"

// Representation of U(1) algebra
u1::u1(FloatType theta){

  field.push_back(theta);
  size=1;

}

u1::u1(std::vector<FloatType> theta){
  if(theta.size()!=1){
    fprintf(stderr,"size is %ld\n", theta.size());    
    throw std::invalid_argument("FundHiggs has the wrong size.");
  } else{
    field=theta;
    size=theta.size();
  }
}




// Representation of SU(2) algebra
su2::su2(std::vector<FloatType> A){
  if(A.size()!=3){
    fprintf(stderr,"size is %ld\n", A.size());    
    throw std::invalid_argument("FundHiggs has the wrong size.");
  } else{
    field=A;
    size=A.size();
  }
}
