// group.cpp

#include <cmath>
#include <complex>
#include <vector>
using namespace std::complex_literals;

#include "precisions.h"
#include "scalar.h"
#include "group.h"
#include "vector.h"



// Action on complex number
std::complex<FloatType> U1::action(const std::complex<FloatType>& a){
  return a*(cos(alg[0])+1.0i*sin(alg[0]));
}

U1 U1::inverse(){
  return U1(-alg);
}

U1 U1::groupMultiply(const U1& U){
  return U1(alg+U.toAlgebra());
}




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

  FloatType a = algebraNorm();
  std::complex<FloatType> alpha = cos(a) + 1.0i*sin(a)*alg[2]/a;
  std::complex<FloatType> beta = sin(a)*(alg[1] + 1.0i*alg[0])/a;

  // Adjoint action of the SU(2) group on 3-d Pauli vector
  std::vector<FloatType> hprime_array;
  hprime_array.push_back(std::real(alpha*alpha-beta*beta)*h[0] + std::imag(alpha*alpha + beta*beta)*h[1] - 2.0*std::real(alpha*beta)*h[2]);
  hprime_array.push_back(std::imag(beta*beta-alpha*alpha)*h[0] + std::real(alpha*alpha + beta*beta)*h[1] + 2.0*std::imag(alpha*beta)*h[2]);
  hprime_array.push_back(2.0*std::real(alpha*std::conj(beta))*h[0] + 2.0*imag(alpha*std::conj(beta))*h[1] + (std::norm(alpha) - std::norm(beta))*h[2]);
  
  return AdHiggs(hprime_array);
  
}

// Group element corresponding to the cross product of su(2) algebra elements taken as Pauli vectors
SU2 SU2::crossProduct(const SU2& U){

  FloatType a = algebraNorm();
  FloatType a2 = U.algebraNorm();
  su2 alg2 = U.toAlgebra();

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

// SU(2) group multiplication expressed from its effect on su(2) elements
SU2 SU2::groupMultiply(const SU2& U){

  FloatType a = algebraNorm();
  FloatType a2 = U.algebraNorm();
  su2 alg2 = U.toAlgebra();

  // Group multiplication from operations on algebra
  FloatType cosc = cos(a)*cos(a2) - (alg*alg2)*sin(a)*sin(a2)/(a*a2);
  su2 s = alg*sin(a)*cos(a2)/a + alg2*cos(a)*sin(a2)/a2 - (crossProduct(U).alg)*sin(a)*sin(a2)/(a*a2);
  FloatType snorm = s.norm();

  if(snorm!=0) s = s/snorm;
  
  return SU2(s*atan2(snorm,cosc));

}

// Print the matrix representation of SU(2) elements
void SU2::matrixForm(){

  FloatType a = algebraNorm();
  std::complex<FloatType> alpha = cos(a) + 1.0i*sin(a)*alg[2]/a;
  std::complex<FloatType> beta = sin(a)*(alg[1] + 1.0i*alg[0])/a;
  printf("\n\n %f + i %f \t %f + i %f \n", std::real(alpha), std::imag(alpha), std::real(beta), std::imag(beta));
  printf(" %f + i %f \t %f + i %f \n\n", -std::real(beta), std::imag(beta), std::real(alpha), -std::imag(alpha));

}

// Print the matrix representation of SU(2) elements, for comparison with groupMultiply
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
