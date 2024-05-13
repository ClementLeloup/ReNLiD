// timeslice.h

#ifndef _TIMESLICE_H_
#define _TIMESLICE_H_

#include "precisions.h"
#include "lattice.h"


// Store lattice information at constant time
class TimeSlice : public VertexLattice
{
public:
  TimeSlice();
  TimeSlice(const TimeSlice& ts);
  ~TimeSlice();
  // Vertex& operator()(const VertexIndex& ind); // get the vertex at index ind from the lattice
  // const Vertex& operator()(const VertexIndex& ind) const;
  TypeFLoat Energy();
  TypeFloat Gauss();
  TypeFloat* PowerSpectrum();
  void WhiteNoise();
  void init(TypeFloat* powerspectrum);
  void init(TimeSlice ts);
  void save();
  
  
protected:
  void alloc();
  void free();
  
  Vertex* vertices;
  TypeFloat a; // Scale factor
};



#endif
