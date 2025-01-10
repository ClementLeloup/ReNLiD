// evolver.h

#ifndef _EVOLVER_H_
#define _EVOLVER_H_

#include "precisions.h"
#include "lattice.h"
#include "timeslice.h"


// Initialize the lattice, evolve it in time, performs sanity checks and save data
class Evolver
{
public:
  Evolver();
  uint getStep();
  FloatType getEnergy();
  FloatType getInflaton();
  FloatType getdInflaton();
  void setEnergy(const FloatType E);
  void setInflaton(const FloatType phi);
  void setdInflaton(const FloatType dphi);
  void setAverage(const FloatType E, const FloatType phi, const FloatType dphi);
  void setAverage(TimeSlice ts);  
  TimeSlice LatticeToTimeSlice();
  void TimeSliceToLattice(const TimeSlice& ts, const uint s);
  void evolveSlice(uint i1);
  void evolveLattice();
  int goToNext();
  void save();

  // MPI communication processes
  // void commRight(VertexLattice &v, int rankTo, int rankFrom);

  
protected:
  uint step;
  FloatType inflaton; // vev of the inflaton field
  FloatType dinflaton; // vev of the derivative wrt time of the inflaton field
  FloatType energy;
  
  VertexLattice v1;
  VertexLattice v2;
  VertexLattice v3;
};


#endif
