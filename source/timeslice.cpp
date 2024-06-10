// timeslice.cpp

#include "precisions.h"
#include "lattice.h"
#include "timeslice.h"
#include <random>
#include <complex.h>
#include <fftw3.h>
#include <cmath>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <functional>

// Store lattice information at constant time
TimeSlice::TimeSlice(){  
  alloc();
  scaleFactor = 1.0;
}

TimeSlice::TimeSlice(const TimeSlice& ts){
  alloc();
  std::memcpy(vertices, ts.vertices, N*N*N*sizeof(Vertex));
  std::memcpy(dvertices, ts.dvertices, N*N*N*sizeof(Vertex));
}

TimeSlice::~TimeSlice(){
  free();
}

FloatType TimeSlice::Energy(){

  // Need implementation
  return 1.0;
  
}

FloatType TimeSlice::Gauss(){

  // Need implementation
  return 1.0;
  
}

// Now with complex fourier transform
std::vector<FloatType> TimeSlice::PowerSpectrum(const std::vector<FloatType> k){
  
  fftw_complex *in, *out;
  fftw_plan p;
  std::vector<FloatType> Pk(k.size(), 0.0);
  std::vector<uint> density(k.size(), 0);

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // input real space data
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // output Fourier space data
  p = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE); // FFT plan

  // std::memcpy(in, vertices, N*N*N*sizeof(fftw_complex)); // copying the data from the lattice to in
  // printf("trucmuche\n");

  for(int i=0; i<N*N*N;i++){
    in[i][0] = vertices[i].getRealScalar()[0];
  }
  
  fftw_execute(p); // FFT

  /* Compute PowerSpectrum */
  for(int i1=0; i1<N; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){

	int index = i3 + i2*N + i1*N*N; // Flatten the indices
        FloatType k1 = (i1 <= N/2) ? i1 : i1-N;
	FloatType k2 = (i2 <= N/2) ? i2 : i2-N;
	FloatType k3 = (i3 <= N/2) ? i3 : i3-N;
	// FloatType ki = std::sqrt(i1*i1 + i2*i2 + i3*i3);  // Get multipole
	FloatType ki = 2.0*M_PI*std::sqrt(k1*k1 + k2*k2 + k3*k3)/(scaleFactor*N);  // Get multipole
	
	uint kindex = 0;
	
	// Search for interval where ki falls in k array
	// Can be optimized
	if(k[0] <= ki && ki <= k[k.size()-1]){
	  while(k[kindex] <= ki && kindex < k.size()){ // k_size TBD
	    kindex++;
	  }
	  density[kindex-1] += 1; // increment density count
	  Pk[kindex-1] += (out[index][0]*out[index][0] + out[index][1]*out[index][1])/(N*N*N); // Check if no square root
	}
	
      }
    }
  }

  // Average over |k| shell
  for(int kindex=0; kindex < k.size(); kindex++){
    if(density[kindex] != 0){Pk[kindex] = Pk[kindex]/density[kindex];}
  }
  
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);

  return Pk;
  
}

void TimeSlice::WhiteNoise(){
  
  std::default_random_engine generator;
  std::normal_distribution<FloatType> distribution(0.0,1.0);

  for (int i=0; i<N*N*N; ++i) {
    FloatType number = distribution(generator);
    vertices[i] = RealScalar(std::vector<FloatType>(1,number));
    dvertices[i] = RealScalar(std::vector<FloatType>(1,number));
  }
  
}

// Be careful, this operates in Fourier space
fftw_complex* TimeSlice::aOperator(){

  fftw_complex *ak;
  ak = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // random realization for the stochastic analog of the annihilation operator
  
  std::default_random_engine generator;
  std::uniform_real_distribution<FloatType> distribution(0.0,1.0);

  for (int k=0; k<N*N*N; ++k) {
    FloatType Xk = distribution(generator);
    FloatType Yk = distribution(generator);
    FloatType ak_mod = std::sqrt(-0.5*log(Xk));
    ak[k][0] = ak_mod*cos(2.0*M_PI*Yk);
    ak[k][1] = ak_mod*sin(2.0*M_PI*Yk);
  }

  return ak;
  
}

void TimeSlice::init(std::function<FloatType(FloatType)> sqrtPk){

  WhiteNoise();

  fftw_complex *in, *out;
  fftw_plan p_fwd, p_bwd;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // input real space data
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // output Fourier space data
  p_fwd = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE); // FFT plan, can maybe reuse a previous one
  p_bwd = fftw_plan_dft_3d(N, N, N, out, in, FFTW_BACKWARD, FFTW_ESTIMATE); // FFT plan, can maybe reuse a previous one

  for(int i=0; i<N*N*N;i++){
    in[i][0] = vertices[i].getRealScalar()[0];
  }
  
  fftw_execute(p_fwd); // FFT

  // Multiply by Pk^1/2
  for(int i1=0; i1<N; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){

	int index = i3 + i2*N + i1*N*N; // Flatten the indices
	FloatType k1 = (i1 <= N/2) ? i1 : i1-N;
	FloatType k2 = (i2 <= N/2) ? i2 : i2-N;
	FloatType k3 = (i3 <= N/2) ? i3 : i3-N;
	// FloatType ki = std::sqrt(i1*i1 + i2*i2 + i3*i3);  // Get multipole
	FloatType ki = 2.0*M_PI*std::sqrt(k1*k1 + k2*k2 + k3*k3)/N;  // Get multipole
	
	out[index][0] *= sqrtPk(ki);
	out[index][1] *= sqrtPk(ki);
      }
    }
  }

  fftw_execute(p_bwd); // inverse FFT  

  // Update the time slice
  for(int i=0; i<N*N*N;i++){
    std::vector<FloatType> psi(1, in[i][0]/(N*N*N));
    vertices[i].setRealScalar(RealScalar(psi));
  }
  
  fftw_destroy_plan(p_fwd);
  fftw_destroy_plan(p_bwd);  
  fftw_free(in);
  fftw_free(out);
  
}

// Initial condition using stochastic equivalent of creation and annihilation operators
void TimeSlice::init(fftw_complex* ak, std::function<std::complex<FloatType>(FloatType)> uk, std::function<std::complex<FloatType>(FloatType)> duk){
  
  fftw_complex *in, *out;
  fftw_plan p_bwd;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // output Fourier space data  
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // output Fourier space data
  din = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // output Fourier space time derivative data  
  dout = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // output Fourier space time derivative data
  p_bwd = fftw_plan_dft_3d(N, N, N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); // FFT plan, can maybe reuse a previous one
  dp_bwd = fftw_plan_dft_3d(N, N, N, din, dout, FFTW_BACKWARD, FFTW_ESTIMATE); // FFT plan, can maybe reuse a previous one

  // Multiply by Pk^1/2
  for(int i1=0; i1<N; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){
	
	int index = i3 + i2*N + i1*N*N; // Flatten the indices
	FloatType k1 = (i1 <= N/2) ? i1 : i1-N;
	FloatType k2 = (i2 <= N/2) ? i2 : i2-N;
	FloatType k3 = (i3 <= N/2) ? i3 : i3-N;
	// FloatType ki = std::sqrt(i1*i1 + i2*i2 + i3*i3);  // Get multipole
	FloatType ki = 2.0*M_PI*std::sqrt(k1*k1 + k2*k2 + k3*k3)/N;  // Get multipole
	int opp_i1 = (i1 == 0) ? 0 : N - i1;
	int opp_i2 = (i2 == 0) ? 0 : N - i2;
	int opp_i3 = (i3 == 0) ? 0 : N - i3;
	int opp_index = opp_i3 + opp_i2*N + opp_i1*N*N; // index of -k, needed because the scalar field is real

	// Spectral decomposition of the scalar field
	in[index][0] = ak[index][0]*std::real(uk(ki)) - ak[index][1]*std::imag(uk(ki)) + ak[opp_index][0]*std::real(uk(ki)) - ak[opp_index][1]*std::imag(uk(ki));
	in[index][1] = ak[index][1]*std::real(uk(ki)) + ak[index][0]*std::imag(uk(ki)) - ak[opp_index][1]*std::real(uk(ki)) - ak[opp_index][0]*std::imag(uk(ki));
	din[index][0] = ak[index][0]*std::real(duk(ki)) - ak[index][1]*std::imag(duk(ki)) + ak[opp_index][0]*std::real(duk(ki)) - ak[opp_index][1]*std::imag(duk(ki));
	din[index][1] = ak[index][1]*std::real(duk(ki)) + ak[index][0]*std::imag(duk(ki)) - ak[opp_index][1]*std::real(duk(ki)) - ak[opp_index][0]*std::imag(duk(ki));
      }
    }
  }
  
  fftw_execute(p_bwd); // inverse FFT
  fftw_execute(dp_bwd); // inverse FFT

  // Update the time slice
  for(int i=0; i<N*N*N;i++){
    std::vector<FloatType> psi(1, out[i][0]/std::sqrt(N*N*N));
    std::vector<FloatType> dpsi(1, dout[i][0]/std::sqrt(N*N*N));
    vertices[i].setRealScalar(RealScalar(psi));
    dvertices[i].setRealScalar(RealScalar(dpsi));
  }
  
  fftw_destroy_plan(p_bwd);  
  fftw_destroy_plan(dp_bwd);  
  fftw_free(in);
  fftw_free(out);
  fftw_free(din);
  fftw_free(dout);
  
}


void TimeSlice::save(std::ofstream& outstream){
  
  if(outstream.is_open()){

    for(int i1=0; i1<N; i1++){
      for(int i2=0; i2<N; i2++){
	for(int i3=0; i3<N; i3++){
	  outstream << i1 << "\t" << i2 << "\t" << i3 << "\t" << vertices[i3 + i2*N + i1*N*N].getRealScalar()[0] << "\t" << dvertices[i3 + i2*N + i1*N*N].getRealScalar()[0] << "\n";
	}
      }
    }
  } else{
    throw std::invalid_argument("Output file is not open.");
  }
}



// protected functions
void TimeSlice::alloc(){  
  vertices = new Vertex[N*N*N];
  dvertices = new Vertex[N*N*N];
}

void TimeSlice::free(){
  delete [] vertices;
  delete [] dvertices;
  vertices = NULL;
  dvertices = NULL;
}
