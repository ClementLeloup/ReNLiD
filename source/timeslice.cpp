// timeslice.cpp

#include "precisions.h"
#include "lattice.h"
#include "timeslice.h"
#include <random>
#include <complex.h>
#include <fftw3.h>


// Store lattice information at constant time
TimeSlice::TimeSlice(){
  alloc();
}

TimeSlice::TimeSlice(const TimeSlice& ts){
  alloc();
  std::memcpy(vertices, ts.vertices, N*N*N*sizeof(Vertex));
}

TimeSlice::~TimeSlice(){
  free();
}

TypeFloat TimeSlice::Energy(){

}

// Now with complex fourier transform
TypeFloat* TimeSlice::PowerSpectrum(){

  fftw_complex *in, *out;
  fftw_plan p;
  TypeFloat* k;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N*N*N);
  p = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  memcpy(in, vertices, N*N*N*sizeof(fftw_complex));
  
  fftw_execute(p);
  for(int i1=0; i<N, i++){
    for(int i2=0; i<N, j++){
      for(int i3=0; i<N, k++){

	
	
      }
    }
  }

  
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
  
}

TypeFloat TimeSlice::Energy(){

}


void TimeSlice::Whitenoise(){
  
  std::default_random_engine generator;
  std::normal_distribution<FloatType> distribution(0.0,1.0);

  for (int i=0; i<N*N*N; ++i) {
    FloatType number = distribution(generator);
    vertices[i] = number;
  }

}

void TimeSlice::init(TypeFloat* powerspectrum){

  a = 1.0;
  this.Whitenoise();
  
  
}

// protected functions
void TimeSlice::alloc(){
  vertices = new Vertex[N*N*N]; // size is the number of MPI processes
}

void TimeSlice::free(){
  delete [] vertices;
  vertices = NULL;
}
