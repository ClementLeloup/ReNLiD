// lattice. cpp

#include <cstring>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
namespace mpi = boost::mpi;

#include "precisions.h"
#include "lattice.h"



// Vertices
// To be extended to multiple or multi-dimensional fields
RealScalar Vertex::getRealScalar(){
  return phi;
}

void Vertex::setRealScalar(RealScalar psi){
  phi = psi;
}

Vertex Vertex::operator+(Vertex v){
  Vertex vsum;
  vsum.setRealScalar(phi+v.getRealScalar());
  return vsum;
}


// Indices of lattice vertices
VertexIndex::VertexIndex(uint i, uint j, uint k){
    ind[0] = i % N; ind[1] = j % N; ind[2] = k % N;
}

VertexIndex::VertexIndex(const VertexIndex& vi){
    ind[0] = vi.ind[0]; ind[1] = vi.ind[1]; ind[2] = vi.ind[2];
}

uint VertexIndex::getIndx() const {
  return ind[0];
}

uint VertexIndex::getIndyz() const {
  return N*ind[1]+ind[2];
}

VertexIndex VertexIndex::operator+(Axis a) const {
  VertexIndex adjacent = (*this);
  adjacent.ind[a] = (ind[a] + 1) % N;
  return adjacent;
}

VertexIndex VertexIndex::operator-(Axis a) const {
  VertexIndex adjacent = (*this);
  adjacent.ind[a] = (ind[a] + N - 1) % N; // Add N to ensure this is positive
  return adjacent;
}


// Indices of lattice edges
uint EdgeIndex::getSlice() const {
  return ind[0]; // Maybe something to modify here
}

uint EdgeIndex::getDir() const {
  return dir;
}


// Vertices of the lattice
// public functions
VertexLattice::VertexLattice(){
  alloc();
}

VertexLattice::VertexLattice(const VertexLattice& vl){
  alloc();
  std::memcpy(vertices, vl.vertices, (N/size + 2)*N*N*sizeof(Vertex));
}

VertexLattice::VertexLattice(const Vertex* v){
  alloc();
  // std::memcpy(vertices, v, (N/size + 2)*N*N*sizeof(Vertex));
  std::copy(v, v+(N/size + 2)*N*N, vertices); // memcpy does not work for some reason...
}

VertexLattice::~VertexLattice(){
  free();
}

Vertex& VertexLattice::operator()(const VertexIndex& ind) {
  return getSlice(ind.getIndx())[ind.getIndyz()];
}

const Vertex& VertexLattice::operator()(const VertexIndex& ind) const {
  return getSlice(ind.getIndx())[ind.getIndyz()];
}

VertexLattice& VertexLattice::operator=(const VertexLattice& vl){
  free();
  alloc();
  std::memcpy(vertices, vl.vertices, (N/size + 2)*N*N*sizeof(Vertex));
  return *this;
}

void VertexLattice::init(){

}

// MPI communication processes
void VertexLattice::send(int rankTo){
  
  // Send last x slice
  if (rankTo == (rank+1) % size){
    requests[2] = world.isend(rankTo, 2*rank, getSlice((rank+1)*N/size - 1), N*N);
  }
  // Send first x slice
  if (rankTo == (size+rank-1) % size){
    requests[0] = world.isend(rankTo, 2*rank+1, getSlice(rank*N/size), N*N);
  }
}

void VertexLattice::receive(int rankFrom){

  // Receive last x slice of previous portion
  if (rankFrom == (size+rank-1) % size){
    requests[1] = world.irecv(rankFrom, 2*rankFrom, getSlice((rankFrom+1)*N/size-1), N*N);
  }
  // Receive first x slice of next portion
  if (rankFrom == (rank+1) % size){
    requests[3] = world.irecv(rankFrom, 2*rankFrom+1, getSlice(rankFrom*N/size), N*N);
  }
}

void VertexLattice::wait(){

  // Wait until all MPI communication is done
  mpi::wait_all(requests, requests+4);

}


// protected functions
void VertexLattice::alloc(){
  vertices = new Vertex[(N/size + 2)*N*N]; // size is the number of MPI processes
}

void VertexLattice::free(){
  delete [] vertices;
  // vertices = NULL;
}

Vertex* VertexLattice::getSlice(uint ind0){
  return vertices + N*N*((ind0 + N + 1 - rank*N/size) % N);
}

const Vertex* VertexLattice::getSlice(uint ind0) const {
  return vertices + N*N*((ind0 + N + 1 - rank*N/size) % N);
}





// Need to implement EdgeLattice functions, analogous to VertexLattice
