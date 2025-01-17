// timeslice.cpp

#include <random>
#include <complex.h>
#include <fftw3.h>
#include <cmath>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <functional>

#include "precisions.h"
#include "lattice.h"
#include "timeslice.h"
#include "scalar.h"



// Store lattice information at constant time
TimeSlice::TimeSlice(){  
  alloc();
  scaleFactor = 1.0;
}

TimeSlice::TimeSlice(const TimeSlice& ts){
  scaleFactor = 1.0;
  alloc();
  std::memcpy(vertices, ts.vertices, N*N*N*sizeof(Vertex));
  std::memcpy(dvertices, ts.dvertices, N*N*N*sizeof(Vertex));
}

TimeSlice::TimeSlice(const Vertex* v, const Vertex* dv){
  scaleFactor = 1.0;
  alloc();
  // std::memcpy(vertices, v, N*N*N*sizeof(Vertex));
  // std::memcpy(dvertices, dv, N*N*N*sizeof(Vertex));
  std::copy(v, v+N*N*N, vertices);
  std::copy(dv, dv+N*N*N, dvertices);
}


TimeSlice::~TimeSlice(){
  free();
}

Vertex* TimeSlice::getField(){
  return vertices;
}

Vertex* TimeSlice::getdField(){
  return dvertices;
}

void TimeSlice::setField(const Vertex* v){
  std::copy(v, v+N*N*N, vertices);
}

void TimeSlice::setdField(const Vertex* dv){
  std::copy(dv, dv+N*N*N, dvertices);
}

Vertex TimeSlice::getFieldAtInd(const VertexIndex& ind) const {
  return vertices[ind.getIndyz() + N*N*ind.getIndx()];
}

Vertex TimeSlice::getdFieldAtInd(const VertexIndex& ind) const {
  return dvertices[ind.getIndyz() + N*N*ind.getIndx()];
}

// Compute the vacuum expectation value of fields in the timeslice
FloatType TimeSlice::vev() const {

  FloatType phi = 0.0;
  
  for(int i1=0; i1<N; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){
	VertexIndex ind(i1, i2, i3); 
	phi += getFieldAtInd(ind).getRealScalar()[0];
      }
    }
  }

  return phi/(N*N*N);
  
}

// Compute the vacuum expectation value of fields derivatives in the timeslice
FloatType TimeSlice::dvev() const {

  FloatType dphi = 0.0;
  
  for(int i1=0; i1<N; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){
	VertexIndex ind(i1, i2, i3); 
	dphi += getdFieldAtInd(ind).getRealScalar()[0];
      }
    }
  }

  return dphi/(N*N*N);
  
}

// Compute total energy of the timeslice
FloatType TimeSlice::Energy() const {

  FloatType E = 0.0;
  
  for(int i1=0; i1<N; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){
	VertexIndex ind(i1, i2, i3); 
	RealScalar phi = getFieldAtInd(ind).getRealScalar();
	RealScalar dphidt = getdFieldAtInd(ind).getRealScalar();
	RealScalar dphidx = (phi - getFieldAtInd(ind-x).getRealScalar())/dx; // Taking left or right derivatives give similar results
	RealScalar dphidy = (phi - getFieldAtInd(ind-y).getRealScalar())/dx;
	RealScalar dphidz = (phi - getFieldAtInd(ind-z).getRealScalar())/dx;
	E += 0.5*(dphidt*dphidt + dphidx*dphidx + dphidy*dphidy + dphidz*dphidz) + phi.V();
      }
    }
  }

  return E;
  
}

// Compute total energy, vev and dvev of the timeslice
std::vector<FloatType> TimeSlice::Average() {

  std::vector<FloatType> avg(3, 0.0);
  for(int i1=0; i1<N; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){
	VertexIndex ind(i1, i2, i3); 
	RealScalar phi = getFieldAtInd(ind).getRealScalar();
	RealScalar dphidt = getdFieldAtInd(ind).getRealScalar();
	RealScalar dphidx = (phi - getFieldAtInd(ind-x).getRealScalar())/dx; // Taking left or right derivatives give similar results
	RealScalar dphidy = (phi - getFieldAtInd(ind-y).getRealScalar())/dx;
	RealScalar dphidz = (phi - getFieldAtInd(ind-z).getRealScalar())/dx;
	avg[0] += 0.5*(dphidt*dphidt + dphidx*dphidx + dphidy*dphidy + dphidz*dphidz) + phi.V(); // energy
	avg[1] += phi[0]; // inflaton
	avg[2] += dphidt[0]; // dinflaton
      }
    }
  }

  avg[1] /= (N*N*N);
  avg[2] /= (N*N*N);
  
  return avg;
  
}


FloatType TimeSlice::Gauss() const {

  // Need implementation when including gauge fields
  
  return 0.0;
  
}

// Calculate powerspectrum of the field spatial distribution
std::vector<FloatType> TimeSlice::PowerSpectrum(const std::vector<FloatType> k){

  FloatType phi0 = vev();
  fftw_complex *in, *out;
  fftw_plan p;
  std::vector<FloatType> Pk(k.size(), 0.0);
  std::vector<int> density(k.size(), 0);

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // input real space data
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // output Fourier space data
  p = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE); // FFT plan
  
  // initialize input and output arrays for Fourier transform
  for(int i=0; i<N*N*N;i++){
    in[i][0] = vertices[i].getRealScalar()[0] - phi0;
    in[i][1] = 0.;
    out[i][0] = 0.;
    out[i][1] = 0.;    
  }
  
  fftw_execute(p); // FFT

  /* Compute PowerSpectrum */
  for(int i1=0; i1<N; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){

	uint index = i3 + i2*N + i1*N*N; // Flatten the indices
        FloatType k1 = (i1 <= N/2) ? i1 : i1-N;
	FloatType k2 = (i2 <= N/2) ? i2 : i2-N;
	FloatType k3 = (i3 <= N/2) ? i3 : i3-N;
	FloatType ki = 2.0*M_PI*std::sqrt(k1*k1 + k2*k2 + k3*k3)/(dx*scaleFactor*N);  // Get multipole
	
	uint kindex = 0;
	
	// Search for interval where ki falls in k array
	// Can be optimized
	if(k[0] <= ki && ki <= k[k.size()-1]){
	  while(kindex < k.size()-1 && k[kindex] + k[kindex+1] < 2.0*ki){ // k_size TBD
	    kindex++;
	  }
	  density[kindex] += 1; // increment density count
	  Pk[kindex] += (out[index][0]*out[index][0] + out[index][1]*out[index][1])/(N*N*N);
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

// Populate the time slice with a white noise realization
void TimeSlice::WhiteNoise(){
  
  std::default_random_engine generator;
  std::normal_distribution<FloatType> distribution(0.0,1.0);

  for (int i=0; i<N*N*N; ++i) {
    FloatType number = distribution(generator);
    vertices[i].setRealScalar(RealScalar(std::vector<FloatType>(1,number)));
    dvertices[i].setRealScalar(RealScalar(std::vector<FloatType>(1,number)));
  }
  
}

// Random creation and annihilation fields, easier to generate the right initial conditions for the fields and their derivatives
fftw_complex* TimeSlice::aOperator(){

  // Be careful, this operates in Fourier space
  fftw_complex *ak;
  ak = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N); // random realization for the stochastic analog of the annihilation operator
  
  std::default_random_engine generator;
  std::uniform_real_distribution<FloatType> distribution(0.0,1.0);

  // Generate a random realization of creation and annihilation fields
  for (int k=0; k<N*N*N; ++k) {
    FloatType Xk = distribution(generator);
    FloatType Yk = distribution(generator);
    FloatType ak_mod = std::sqrt(-0.5*log(Xk));
    ak[k][0] = ak_mod*cos(2.0*M_PI*Yk);
    ak[k][1] = ak_mod*sin(2.0*M_PI*Yk);
  }

  return ak;
  
}

// Initialize the initial condition to a random realization with correct statistical properties, from the field powerspectrum square root
void TimeSlice::init(RealScalar phi, std::function<FloatType(FloatType)> sqrtPk){

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
	FloatType ki = 2.0*M_PI*std::sqrt(k1*k1 + k2*k2 + k3*k3)/(dx*scaleFactor*N);  // Get multipole
	
	out[index][0] *= sqrtPk(ki);
	out[index][1] *= sqrtPk(ki);
      }
    }
  }

  fftw_execute(p_bwd); // inverse FFT  

  // Update the time slice
  for(int i=0; i<N*N*N;i++){
    std::vector<FloatType> psi(1, in[i][0]/(N*N*N));
    vertices[i].setRealScalar(phi+RealScalar(psi));
  }
  
  fftw_destroy_plan(p_fwd);
  fftw_destroy_plan(p_bwd);  
  fftw_free(in);
  fftw_free(out);
  
}

// Initial condition using stochastic equivalent of creation and annihilation operators
void TimeSlice::init(RealScalar phi, RealScalar dphi, fftw_complex* ak, std::function<std::complex<FloatType>(FloatType)> uk, std::function<std::complex<FloatType>(FloatType)> duk){
  
  fftw_complex *in, *out;
  fftw_complex *din, *dout;  
  fftw_plan p_bwd, dp_bwd;

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
	FloatType ki = 2.0*M_PI*std::sqrt(k1*k1 + k2*k2 + k3*k3)/(dx*scaleFactor*N);  // Get multipole
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
    vertices[i].setRealScalar(phi+RealScalar(psi));
    dvertices[i].setRealScalar(dphi+RealScalar(dpsi));
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
  // vertices = NULL;
  // dvertices = NULL;
}
