// ReNLiD.cpp

#include <stdio.h>
#include <iostream>
#include <typeinfo>

// Import complex type
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
using namespace std::complex_literals;

#include "precisions.h"
#include "field.h"
#include "lie.h"
#include "scalar.h"
#include "group.h"




/////////////////////////////////////////////////////////////////////////
/////////////////////////// Testing functions ///////////////////////////
/////////////////////////////////////////////////////////////////////////

void test_FundHiggs(){

  // std::complex<FloatType> test_c[2] = {1.5+1.0i, 3.0+1.0i};
  // std::complex<FloatType> test_c2[2] = {1.5+1.5i, 3.0+1.5i};
  std::vector<std::complex<FloatType>> test_c = {1.5+1.0i, 3.0+1.0i};
  std::vector<std::complex<FloatType>> test_c2 = {1.5+1.5i, 3.0+1.5i};

  
  printf("\nHiggs:\n");
  FundHiggs higgs(test_c);
  higgs.print();
  // higgs = higgs*(2.0+3.0i);
  (higgs*(2.0+3.0i)).print();
  higgs.print();

  printf("\nHiggs2:\n");
  FundHiggs higgs2(test_c2);
  higgs2.print();
  (higgs2/(2.0+1.0i)).print();
  higgs2.print();
  
  (higgs+higgs2).print();
  
  // Scalar product
  std::cout << "Norm squared = " << higgs.conj()*higgs << "\t but also = " << higgs.normSquared() << "\n";
  std::cout << "Scalar product = " << higgs.conj()*higgs2 << "\n";
  
  
};

void test_AdHiggs(){

  // FloatType test_c[3] = {1.5, 3.0, 4.5};
  // FloatType test_c2[3] = {1.0, 2.0, 3.0};
  std::vector<FloatType> test_c = {1.5, 3.0, 4.5};
  std::vector<FloatType> test_c2 = {1.0, 2.0, 3.0};

  printf("\nHiggs:\n");
  AdHiggs higgs(test_c);
  higgs.print();
  printf("2 x higgs =\n");
  //higgs = higgs*2.0;
  (higgs*2.0).print();
  higgs.print();
  
  printf("\nHiggs2:\n");
  AdHiggs higgs2(test_c2);
  higgs2.print();
  printf("higgs2 / 2 =\n");
  (higgs2/2.0).print();
  higgs2.print();

  printf("\nHiggs + Higgs2 =\n");
  (higgs+higgs2).print();
  // printf("\nHiggs x Higgs2 =\n");
  // (higgs*higgs2).print();


  // Scalar product
  std::cout << "Norm squared = " << higgs*higgs << "\t but also = " << higgs.normSquared() << "\n";
  std::cout << "Scalar product = " << higgs*higgs2 << "\n";
  
  
};

void test_U1(){

  // FloatType theta1[1] = {2.5};
  // FloatType theta2[1] = {1.5};
  u1 A1(2.5);
  u1 A2(1.5);
  printf("\nu1: \n");
  A1.print();
  A2.print();
  
  U1 g1(A1);
  U1 g2(A2);
  printf("\nU1: \n");
  g1.print();
  g2.print();

  FloatType a = 1.5;
  printf("\ng1.a = %f + i%f\n", (g1.action(a)).real(), (g1.action(a)).imag());
  printf("\ng2.a = %f + i%f\n", (g2.action(a)).real(), (g2.action(a)).imag());

  (g1.groupMultiply(g2)).print();
  (g2.groupMultiply(g1)).print();  

}

void test_SU2(){

  std::vector<FloatType> a1 = {2.0, 3.0, 4.0};
  std::vector<FloatType> a2 = {1.5, 2.5, 3.5};
  // std::vector<FloatType> a2 = {1., 1.5, 2};
  su2 A1(a1);
  su2 A2(a2);
  printf("\nsu2: \n");
  // A1.print();
  // A2.print();
  
  SU2 g1(A1);
  SU2 g2(A2);
  printf("\nSU2: \n");
  // g1.print();
  // g2.print();

  std::vector<FloatType> a = {0.5, 1.0, 1.5};
  AdHiggs h(a);
  // h.print();
  AdHiggs truc = g1.action(h);
  // truc.print();
  // printf("truc0 = %f, truc1 = %f, truc2 = %f\n", truc[0], truc[1], truc[2]);
  // h.print();
  // (g2.action(h)).print();

  SU2 g1g2 = g1.groupMultiply(g2);
  SU2 g2g1 = g2.groupMultiply(g1);
 
  g1.matrixForm();
  g2.matrixForm();

  g1g2.matrixForm();
  g2g1.matrixForm();

  g1.productMatrixForm(g2);
  g2.productMatrixForm(g1);
  
  SU2 machin = g1.groupMultiply(g1.inverse());
  machin.print();
  //(machin.action(h)).print();


}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



int main(){


  // printf("Real:\n");
  // test_phi();

  
  // printf("\nComplex:\n");
  // test_phic();

  
  // printf("Fundamental Higgs:\n");
  // test_FundHiggs();

  // printf("\nAdjoint Higgs:\n");
  // test_AdHiggs();  

  // printf("U1 rep:\n");
  // test_U1();

  printf("\nSU2 rep:\n");
  test_SU2();

  return 0;
  
}
