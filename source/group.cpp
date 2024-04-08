// group.cpp

// #include <math.h>
#include <cmath>
#include <complex>
#include <vector>
using namespace std::complex_literals;

#include "precisions.h"
#include "scalar.h"
#include "group.h"
#include "vector.h"


// Representation of U(1) algebra
// u1::u1(FloatType theta){
//   FloatType field_theta[1] = {theta};
//   field=field_theta;
//   size=1;
// }

// u1::u1(FloatType theta){

//   field.push_back(theta);
//   size=1;

// }

// u1::u1(std::vector<FloatType> theta){
//   if(theta.size()!=1){
//     fprintf(stderr,"size is %ld\n", theta.size());    
//     throw std::invalid_argument("FundHiggs has the wrong size.");
//   } else{
//     field=theta;
//     size=theta.size();
//   }
// }


// Representation of U(1) group
// Act on real number
// std::complex<FloatType> U1::act(const FloatType& a){
//   return a*(cos(alg)+i*sin(alg));
// }

// Act on complex number
std::complex<FloatType> U1::action(const std::complex<FloatType>& a){
  return a*(cos(alg[0])+1.0i*sin(alg[0]));
}

U1 U1::inverse(){
  return U1(-alg);
}

U1 U1::groupMultiply(const U1& U){
  return U1(alg+U.toAlgebra());
}




// // Representation of SU(2) algebra
// su2::su2(std::vector<FloatType> A){
//   if(A.size()!=3){
//     fprintf(stderr,"size is %ld\n", A.size());    
//     throw std::invalid_argument("FundHiggs has the wrong size.");
//   } else{
//     field=A;
//     size=A.size();
//   }
// }

// Representation of SU(2) group

//         (0 1)          (0 -i)          (1  0)
// sig_1 = (1 0), sig_2 = (i  0), sig_3 = (0 -1)

// // Constructors
// SU2::SU2(Field A){

//   if(A.getSize()!=3){
//     throw std::invalid_argument("su(2) element has the wrong size.");
//   } else{
//     alg=A;
//   }

// }

// Act on SU2 rep
FundHiggs SU2::action(const FundHiggs& h){

  FloatType a = algebraNorm();
  std::vector<std::complex<FloatType>> hprime_array;
  hprime_array.push_back(cos(a)*h[0] + 1.0i*sin(a)*(alg[0]*h[1]-1.0i*alg[1]*h[1]+alg[2]*h[0])/a);
  hprime_array.push_back(cos(a)*h[1] + 1.0i*sin(a)*(alg[0]*h[0]+1.0i*alg[1]*h[0]-alg[2]*h[1])/a);
  return FundHiggs(hprime_array);
  
}

// Act on SO3 rep
AdHiggs SU2::action(const AdHiggs& h){

  // printf("h0 = %f, h1 = %f, h2 = %f\n", h[0], h[1], h[2]);
  // printf("alg0 = %f, alg1 = %f, alg2 = %f\n", alg[0], alg[1], alg[2]);
  
  FloatType a = algebraNorm();
  // printf("Norm = %f\n", a);
  std::complex<FloatType> alpha = cos(a) + 1.0i*sin(a)*alg[2]/a;
  std::complex<FloatType> beta = sin(a)*(alg[1] + 1.0i*alg[0])/a;

  printf("alpha = %f + i%f\n", std::real(alpha), std::imag(alpha));
  printf("beta = %f + i%f\n", std::real(beta), std::imag(beta));

  // printf("      %f\t%f\t%f\n", std::imag(beta*beta-alpha*alpha), std::real(alpha*alpha + beta*beta), 2.0*std::imag(alpha*beta));
  // printf("Phi = %f\t%f\t%f\n", std::real(beta*beta-alpha*alpha), std::imag(alpha*alpha + beta*beta), 2.0*std::real(alpha*beta));
  // printf("      %f\t%f\t%f\n", 2.0*std::real(alpha*std::conj(beta)), - 2.0*imag(alpha*std::conj(beta)), (std::norm(alpha) - std::norm(beta)));

  // FloatType hprime_array[3];
  // hprime_array[0] = std::imag(beta*beta-alpha*alpha)*h[0] + std::real(alpha*alpha + beta*beta)*h[1] + 2.0*std::imag(alpha*beta)*h[2];
  // hprime_array[1] = std::real(beta*beta-alpha*alpha)*h[0] + std::imag(alpha*alpha + beta*beta)*h[1] + 2.0*std::real(alpha*beta)*h[2];
  // hprime_array[2] = 2.0*std::real(alpha*std::conj(beta))*h[0] - 2.0*imag(alpha*std::conj(beta))*h[1] + (std::norm(alpha) - std::norm(beta))*h[2];

  // std::vector<FloatType> hprime_array;
  // hprime_array.push_back(std::imag(beta*beta-alpha*alpha)*h[0] + std::real(alpha*alpha + beta*beta)*h[1] + 2.0*std::imag(alpha*beta)*h[2]);
  // hprime_array.push_back(std::real(beta*beta-alpha*alpha)*h[0] + std::imag(alpha*alpha + beta*beta)*h[1] + 2.0*std::real(alpha*beta)*h[2]);
  // hprime_array.push_back(2.0*std::real(alpha*std::conj(beta))*h[0] - 2.0*imag(alpha*std::conj(beta))*h[1] + (std::norm(alpha) - std::norm(beta))*h[2]);

  std::vector<FloatType> hprime_array;
  hprime_array.push_back(std::real(alpha*alpha-beta*beta)*h[0] + std::imag(alpha*alpha + beta*beta)*h[1] - 2.0*std::real(alpha*beta)*h[2]);
  hprime_array.push_back(std::imag(beta*beta-alpha*alpha)*h[0] + std::real(alpha*alpha + beta*beta)*h[1] + 2.0*std::imag(alpha*beta)*h[2]);
  hprime_array.push_back(2.0*std::real(alpha*std::conj(beta))*h[0] + 2.0*imag(alpha*std::conj(beta))*h[1] + (std::norm(alpha) - std::norm(beta))*h[2]);
  
  printf("h0 = %f, h1 = %f, h2 = %f\n", hprime_array[0], hprime_array[1], hprime_array[2]);

  return AdHiggs(hprime_array);
  // return hprime_array;
  
}


SU2 SU2::crossProduct(const SU2& U){

  FloatType a = algebraNorm();
  FloatType a2 = U.algebraNorm();

  Field alg2 = U.toAlgebra();

  std::vector<FloatType> algprime_array;
  algprime_array.push_back((alg[1]*alg2[2] - alg[2]*alg2[1]));
  algprime_array.push_back((alg[2]*alg2[0] - alg[0]*alg2[2]));
  algprime_array.push_back((alg[0]*alg2[1] - alg[1]*alg2[0]));

  su2 a1Crossa2(algprime_array);
  return SU2(a1Crossa2);
}

SU2 SU2::inverse(){
  return SU2(-alg);
}

SU2 SU2::groupMultiply(const SU2& U){

  FloatType a = algebraNorm();
  FloatType a2 = U.algebraNorm();

  su2 alg2 = U.toAlgebra();
  
  FloatType cosc = cos(a)*cos(a2) - (alg*alg2)*sin(a)*sin(a2)/(a*a2);
  su2 s = alg*sin(a)*cos(a2)/a + alg2*cos(a)*sin(a2)/a2 - (crossProduct(U).alg)*sin(a)*sin(a2)/(a*a2);
  FloatType snorm = s.norm();

  if(snorm!=0) s = s/snorm;
  
  alg.print();
  alg2.print();
  // (alg*alg2).print();
  s.print();
  printf("tan(c) = %f, c = %f, a1.a2 = %f, |a1| = %f, |a2| = %f, cos(c) = %f and |s| = %f\n", snorm/cosc, atan2(snorm,cosc), alg*alg2, a, a2, cosc, snorm);
  
  return SU2(s*atan2(snorm,cosc));///snorm);

}

void SU2::matrixForm(){


  FloatType a = algebraNorm();
  std::complex<FloatType> alpha = cos(a) + 1.0i*sin(a)*alg[2]/a;
  std::complex<FloatType> beta = sin(a)*(alg[1] + 1.0i*alg[0])/a;
  printf("\n\n %f + i %f \t %f + i %f \n", std::real(alpha), std::imag(alpha), std::real(beta), std::imag(beta));
  printf(" %f + i %f \t %f + i %f \n\n", -std::real(beta), std::imag(beta), std::real(alpha), -std::imag(alpha));

}


void SU2::productMatrixForm(const SU2& U){

  FloatType a = algebraNorm();
  std::complex<FloatType> alpha = cos(a) + 1.0i*sin(a)*alg[2]/a;
  std::complex<FloatType> beta = sin(a)*(alg[1] + 1.0i*alg[0])/a;
  
  FloatType a2 = U.algebraNorm();
  su2 alg2 = U.toAlgebra();
  std::complex<FloatType> alpha2 = cos(a2) + 1.0i*sin(a2)*alg2[2]/a2;
  std::complex<FloatType> beta2 = sin(a2)*(alg2[1] + 1.0i*alg2[0])/a2;

  std::complex<FloatType> prod00 = alpha*alpha2 - beta*std::conj(beta2);
  std::complex<FloatType> prod01 = -std::conj(beta)*alpha2 - std::conj(alpha*beta2);
  std::complex<FloatType> prod10 = alpha*beta2 + beta*std::conj(alpha2);
  std::complex<FloatType> prod11 = -std::conj(beta)*beta2 + std::conj(alpha*alpha2);
  
  printf("\n\n %f + i %f \t %f + i %f \n", std::real(prod00), std::imag(prod00), std::real(prod01), std::imag(prod01));
  printf(" %f + i %f \t %f + i %f \n\n", std::real(prod10), std::imag(prod10), std::real(prod11), std::imag(prod11));
  
}
