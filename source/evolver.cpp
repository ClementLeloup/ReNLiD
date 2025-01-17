// evolver.cpp

#include <cstdlib>
#include <string>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
namespace mpi = boost::mpi;

#include "precisions.h"
#include "scalar.h"
#include "lattice.h"
#include "timeslice.h"
#include "evolver.h"



Evolver::Evolver(){

  step = 0;
  energy = 0;
  inflaton = 1; // initial condition for the inflaton field
  dinflaton = 0; // initial condition for the inflaton derivative wrt time
  
}

uint Evolver::getStep(){
  return step;
}

FloatType Evolver::getEnergy(){
  return energy;
}

FloatType Evolver::getInflaton(){
  return inflaton;
}

FloatType Evolver::getdInflaton(){
  return dinflaton;
}

void Evolver::setEnergy(const FloatType E){
  energy = E;
}

void Evolver::setInflaton(const FloatType phi){
  inflaton = phi;
}

void Evolver::setdInflaton(const FloatType dphi){
  dinflaton = dphi;
}

void Evolver::setAverage(const FloatType E, const FloatType phi, const FloatType dphi){
  energy = E;
  inflaton = phi;
  dinflaton = dphi;
}

void Evolver::setAverage(TimeSlice ts){

  std::vector<FloatType> avg = ts.Average();
  energy = avg[0];
  inflaton = avg[1];
  dinflaton = avg[2];

}

// Gather distributed lattice at given time into single time slice
TimeSlice Evolver::LatticeToTimeSlice(){

  Vertex* vertices;
  Vertex* dvertices;

  if(rank == 0){
    vertices = new Vertex[N*N*N];
    dvertices = new Vertex[N*N*N];
  } else {
    vertices = new Vertex[N*N*N/size];
    dvertices = new Vertex[N*N*N/size];
  }

  // Construct field and field derivative arrays to populate the time slice
  // To be extended to other fields
  for(int i1=0; i1<N/size; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){
	VertexIndex ind(i1+rank*N/size,i2,i3);
	vertices[i3+N*i2+N*N*i1].setRealScalar(v2(ind).getRealScalar());
	dvertices[i3+N*i2+N*N*i1].setRealScalar((v2(ind).getRealScalar() - v1(ind).getRealScalar())/dtau);
      }
    }
  }

  // MPI communication here
  if(rank != 0){
    world.send(0, 2*rank, vertices, N*N*N/size);
    world.send(0, 2*rank+1, dvertices, N*N*N/size);
  } else{
    for(int r=1; r<size; r++){
     world.recv(r, 2*r, &(vertices[r*N*N*N/size]), N*N*N/size);
      world.recv(r, 2*r+1, &(dvertices[r*N*N*N/size]), N*N*N/size);
    }
  }


  // Create and populate the time slice
  TimeSlice ts;
  if(rank == 0){
    ts.setField(vertices);
    ts.setdField(dvertices);
  }

  delete [] vertices;
  delete [] dvertices;
  
  return ts;
  
  
}


// Initialize distributed lattice from time slice
void Evolver::TimeSliceToLattice(const TimeSlice& ts, const uint s){

  step = s;
  energy = ts.Energy();
  inflaton = ts.vev();
  dinflaton = ts.dvev();

  // Construct lattice state at two times from fields and field derivatives in time slice
  // To be extended to other fields
  for(int i1=rank*N/size-1; i1<(rank+1)*N/size+1; i1++){
    for(int i2=0; i2<N; i2++){
      for(int i3=0; i3<N; i3++){
	int i1_pos = (N + i1) % N;//(i1 >= 0) ? i1 : N-i1; // For the case when rank=0 and i1 starts at -1
	VertexIndex ind(i1_pos,i2,i3);
	v1(ind).setRealScalar(ts.getFieldAtInd(ind).getRealScalar());
	v2(ind).setRealScalar(ts.getFieldAtInd(ind).getRealScalar() + ts.getdFieldAtInd(ind).getRealScalar()*dtau);
      }
    }
  }
  
}

// Where the physics is, determine the state of the lattice at the next time step from the state of the lattice in the two previous steps
void Evolver::evolveSlice(uint i1){

  for(int i2=0; i2<N; i2++){
    for(int i3=0; i3<N; i3++){
      VertexIndex ind(i1,i2,i3);

      RealScalar phi1 = v1(ind).getRealScalar();
      RealScalar phi2 = v2(ind).getRealScalar()*2.0;

      // To be extended to other fields
      // v3(ind).setRealScalar(phi2 - phi1 + ((d2phidx2 + d2phidy2 + d2phidz2)/(dx*dx) - (phi2*0.5).dVdphi())*dtau*dtau);
      v3(ind).setRealScalar(phi2 - phi1 + ((v2(ind+x).getRealScalar() + v2(ind-x).getRealScalar()
					    + v2(ind+y).getRealScalar() + v2(ind-y).getRealScalar()
					    + v2(ind+z).getRealScalar() + v2(ind-z).getRealScalar()
					    - phi2*3.0)/(dx*dx) - (phi2*0.5).dVdphi())*dtau*dtau);
    }
  }

  /// Need to add evolution of the EdgeLattice ///
  
}


// Manage computation and MPI communication needed to go to next time step
void Evolver::evolveLattice(){

  int rankLeft = (size + rank - 1) % size;
  int rankRight = (rank + 1) % size;
  
  evolveSlice(rank*N/size); // Evolve leftmost slice first to send to other MPI process
  v3.send(rankLeft);
  v3.receive(rankRight);
  
  evolveSlice((rank+1)*N/size - 1); // Evolve rightmost slice second to send to other MPI process
  v3.send(rankRight);
  v3.receive(rankLeft);
  
  for(int i1=rank*N/size+1; i1<(rank+1)*N/size-1; i1++){
    evolveSlice(i1);
  }
  v3.wait();
  
  /// Need to evolve the EdgeLattice ///
  
}


// Roll over three slices to evolve in time with minimal storage need
int Evolver::goToNext(){
  
  Vertex* temp = v1.vertices;
  v1.vertices = v2.vertices;
  v2.vertices = v3.vertices;
  v3.vertices = temp;

  /// Need to roll over the EdgeLattices ///


  step++;
  
  if(step % 10 == 0){
    // save();
    TimeSlice ts = LatticeToTimeSlice();
    if(rank == 0){
      setAverage(ts.Energy(), ts.vev(), ts.dvev());
    }
    // Need to communicate the status to other MPI processes
  }

  return EXIT_SUCCESS;
}


void Evolver::save(){

  // Need implementation

}
