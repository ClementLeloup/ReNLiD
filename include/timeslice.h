// timeslice.h

#ifndef _TIMESLICE_H_
#define _TIMESLICE_H_

#include "precisions.h"
#include "lattice.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <functional>
#include <complex>
#include <fftw3.h>

// Store lattice information at constant time
class TimeSlice// : public VertexLattice
{
public:
  TimeSlice();
  TimeSlice(const TimeSlice& ts);
  ~TimeSlice();
  Vertex& operator()(const VertexIndex& ind); // get the vertex at index ind from the lattice
  // const Vertex& operator()(const VertexIndex& ind) const;
  FloatType Energy();
  FloatType Gauss();
  std::vector<FloatType> PowerSpectrum(const std::vector<FloatType> k);
  void WhiteNoise();
  fftw_complex* aOperator(); // This operates in Fourier space
  // void init(const std::vector<FloatType> k, const std::vector<FloatType> sqrtPk);
  void init(std::function<FloatType(FloatType)> sqrtPk);
  void init(fftw_complex* ak, std::function<std::complex<FloatType>(FloatType)> uk, std::function<std::complex<FloatType>(FloatType)> duk);
  void init(TimeSlice ts);
  void save(std::ofstream& outstream);
  
  
protected:
  void alloc();
  void free();
  
  Vertex* vertices;
  Vertex* dvertices; // Time derivative of the fields
  FloatType scaleFactor; // Scale factor
};



#endif
