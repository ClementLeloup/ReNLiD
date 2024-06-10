// lattice. cpp

#include "mpi.h"
#include "precisions.h"
#include "lattice.h"
#include <cstring>

// Vertices
RealScalar Vertex::getRealScalar(){
  return phi;
}

void Vertex::setRealScalar(RealScalar psi){
  phi = psi;
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

void VertexLattice::goToNext(VertexLattice& vl1, VertexLattice& vl2, VertexLattice& vl3){
  Vertex* temp = vl1.vertices;
  vl1.vertices = vl2.vertices;
  vl2.vertices = vl3.vertices;
  vl3.vertices = temp;
}

void VertexLattice::init(){

}

// MPI communication processes
void VertexLattice::send(int rankTo){
  // Send last x slice
  if (rankTo == (rank+1 % size)){
    MPI_Isend(getSlice((rank+1)*N/size - 1), N*N*sizeof(Vertex), MPI_BYTE, rankTo, rank, MPI_COMM_WORLD, &rightSend);
  }
  // Send first x slice
  if (rankTo == (rank-1 % size)){
    MPI_Isend(getSlice(rank*N/size), N*N*sizeof(Vertex), MPI_BYTE, rankTo, rank, MPI_COMM_WORLD, &leftSend);
  }
}

void VertexLattice::receive(int rankFrom){
  // Receive last x slice of previous portion
  if (rankFrom == (rank-1 % size)){
    MPI_Irecv(getSlice((rankFrom+1)*N/size - 1), N*N*sizeof(Vertex), MPI_BYTE, rankFrom, rankFrom, MPI_COMM_WORLD, &rightRec);
  }
  // Receive first x slice of next portion
  if (rankFrom == (rank+1 % size)){
    MPI_Irecv(getSlice(rankFrom*N/size), N*N*sizeof(Vertex), MPI_BYTE, rankFrom, rankFrom, MPI_COMM_WORLD, &leftRec);
  }
}

void VertexLattice::wait(){
  MPI_Wait(&leftSend, MPI_STATUS_IGNORE);
  MPI_Wait(&leftRec, MPI_STATUS_IGNORE);
  MPI_Wait(&rightSend, MPI_STATUS_IGNORE);
  MPI_Wait(&rightRec, MPI_STATUS_IGNORE);
}





// protected functions
void VertexLattice::alloc(){
  vertices = new Vertex[(N/size + 2)*N*N]; // size is the number of MPI processes
}

void VertexLattice::free(){
  delete [] vertices;
  vertices = NULL;
}

Vertex* VertexLattice::getSlice(uint ind0){
  return vertices + N*N*((ind0 + N + 1 - rank*N/size) % N);
}

const Vertex* VertexLattice::getSlice(uint ind0) const {
  return vertices + N*N*((ind0 + N + 1 - rank*N/size) % N);
}
