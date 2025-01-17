// timeslice.h

#ifndef _TIMESLICE_H_
#define _TIMESLICE_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <functional>
#include <complex>
#include <fftw3.h>

#include "precisions.h"
#include "lattice.h"
#include "scalar.h"


// Stores lattice information at a given time
class TimeSlice
{
public:
  TimeSlice();
  TimeSlice(const TimeSlice& ts);
  TimeSlice(const Vertex* v, const Vertex* dv);
  ~TimeSlice();

  // Return value of class members
  Vertex* getField(); // get the vertex array
  Vertex* getdField(); // get the dvertex array
  void setField(const Vertex* v);
  void setdField(const Vertex* dv);
  Vertex getFieldAtInd(const VertexIndex& ind) const; // get the field at index ind from the time slice
  Vertex getdFieldAtInd(const VertexIndex& ind) const; // get the field derivative wrt time at index ind from the time slice
  FloatType getScaleFactor();

  FloatType vev() const;
  FloatType dvev() const;
  FloatType Energy() const;
  std::vector<FloatType> Average();
  FloatType Gauss() const;
  std::vector<FloatType> PowerSpectrum(const std::vector<FloatType> k);
  void WhiteNoise();
  fftw_complex* aOperator(); // This operates in Fourier space
  void init(RealScalar phi, std::function<FloatType(FloatType)> sqrtPk);
  void init(RealScalar phi, RealScalar dphi, fftw_complex* ak, std::function<std::complex<FloatType>(FloatType)> uk, std::function<std::complex<FloatType>(FloatType)> duk);
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
