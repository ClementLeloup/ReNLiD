// test_timeslice.cpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <math.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
#include <fftw3.h>
using namespace std;

#include "precisions.h"
#include "timeslice.h"
#include <complex>

int rank;
int size;

FloatType sqrtPk(FloatType k){
  FloatType sqrtPk = (k != 0.0) ? pow(k, -2.0) : 0.0;
  return sqrtPk;
}

// m = 0 and tau = 1 for now
std::complex<FloatType> uk(FloatType k){
  std::complex<FloatType> uk = (k != 0.0) ? exp(-1i*k)/sqrt(2.0*k) : 0.0;
  return uk;
}

// m = 0 for now
std::complex<FloatType> duk(FloatType k){
  std::complex<FloatType> uk = (k != 0.0) ? -1i*k*exp(-1i*k)/sqrt(2.0*k) : 0.0;
  
  return uk;
}

int main(int argc, char **argv){

  std::vector<FloatType> k(N, 0.0);
  for(int i=0;i<N;i++){
    // k[i] = (2.0*M_PI*i)/N;
    k[i] = M_PI*i/(dx*N);
  }

  /*
  printf("Test TimeSlice:\n");
  
  TimeSlice ts;

  printf("Producing white noise\n");
  ts.WhiteNoise();

  printf("Computing power spectrum\n");
  std::vector<FloatType> Pk = ts.PowerSpectrum(k);
  
  std::ofstream file_slice;
  file_slice.open("data/test_timeslice.txt", ios::out);// | ios::app);
  ts.save(file_slice);
  file_slice.close();

  std::ofstream file_pk;
  file_pk.open("data/test_timeslice_Pk.txt", ios::out);// | ios::app);
  for(int i=0;i<N;i++){
    file_pk << k[i] << "\t" << Pk[i] << "\n";
  }
  file_pk.close();

  printf("Producing initial condition with power spectrum\n");
  TimeSlice ts2;
  ts2.init(&sqrtPk);

  std::vector<FloatType> Pk_init = ts2.PowerSpectrum(k);
  
  std::ofstream file_slice_init;
  file_slice_init.open("data/test_timeslice_init.txt", ios::out);// | ios::app);
  ts2.save(file_slice_init);
  file_slice_init.close();

  std::ofstream file_pk_init;
  file_pk_init.open("data/test_timeslice_Pk_init.txt", ios::out);// | ios::app);
  for(int i=0;i<N;i++){
    file_pk_init << k[i] << "\t" << Pk_init[i] << "\n";
  }
  file_pk_init.close();
  */





  {
    printf("Producing initial condition with creation/annihilation\n");
    TimeSlice ts3;  
    fftw_complex* ak = ts3.aOperator();
    ts3.init(ak, &uk, &duk);
    fftw_free(ak);
    
    std::vector<FloatType> Pk_op = ts3.PowerSpectrum(k);
    
    std::ofstream file_slice_op;
    file_slice_op.open("data/test_timeslice_op.txt", ios::out);// | ios::app);
    ts3.save(file_slice_op);
    file_slice_op.close();

    std::ofstream file_pk_op;
    file_pk_op.open("data/test_timeslice_Pk_op.txt", ios::out);// | ios::app);
    for(int i=0;i<k.size();i++){
      file_pk_op << k[i] << "\t" << Pk_op[i] << "\n";
    }
    file_pk_op.close();
  }
  
  return 0;
  
}
