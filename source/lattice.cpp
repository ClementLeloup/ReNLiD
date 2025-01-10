// lattice. cpp

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
namespace mpi = boost::mpi;
// #include "mpi.h"

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

Vertex Vertex::operator+(Vertex v){
  Vertex vsum;
  vsum.setRealScalar(phi+v.getRealScalar());
  // Vertex vsum(phi+v.getRealScalar());
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
  std::copy(v, v+(N/size + 2)*N*N, vertices);
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
  // printf("Send from %d to %d\n", rank, rankTo);
  
  // Send last x slice
  if (rankTo == (rank+1) % size){
    // Vertex* test = getSlice((rank+1)*N/size - 1);
    // for(int i = 0; i< N*N; i++){
    //   printf("Send to right. Element %d in rank %d is %f\n", i, rank, test[i].getRealScalar()[0]);
    // }
       
    // MPI_Isend(getSlice((rank+1)*N/size - 1), N*N*sizeof(Vertex), MPI_BYTE, rankTo, rank, MPI_COMM_WORLD, &rightSend);
    // rightSend = world.isend(rankTo, rank, getSlice((rank+1)*N/size - 1), N*N);
    // printf("Rank is %d and rankTo is %d\n", rank, rankTo);
    // world.send(rankTo, 2*rank, getSlice((rank+1)*N/size - 1), N*N);
    requests[2] = world.isend(rankTo, 2*rank, getSlice((rank+1)*N/size - 1), N*N);
    // printf("2\n");
  }
  // Send first x slice
  if (rankTo == (size+rank-1) % size){
    // Vertex* test = getSlice(rank*N/size - 1);
    // for(int i = 0; i< N*N; i++){
    //   printf("Send to left. Element %d in rank %d is %f\n", i, rank, test[i].getRealScalar()[0]);
    // }
    // printf("I am here, %f \n");
    
    // MPI_Isend(getSlice(rank*N/size), N*N*sizeof(Vertex), MPI_BYTE, rankTo, rank, MPI_COMM_WORLD, &leftSend);
    // leftSend = world.isend(rankTo, rank, getSlice(rank*N/size), N*N);
    // world.send(rankTo, 2*rank+1, getSlice(rank*N/size), N*N);
    requests[0] = world.isend(rankTo, 2*rank+1, getSlice(rank*N/size), N*N);
    // printf("0\n");
  }
}

void VertexLattice::receive(int rankFrom){
  // printf("Receive from %d to %d\n", rankFrom, rank);

  // Receive last x slice of previous portion
  if (rankFrom == (size+rank-1) % size){
    // MPI_Irecv(getSlice((rankFrom+1)*N/size - 1), N*N*sizeof(Vertex), MPI_BYTE, rankFrom, rankFrom, MPI_COMM_WORLD, &rightRec);
    // leftRec = world.irecv(rankFrom, rankFrom, getSlice(rankFrom*N/size-1), N*N);
    // world.recv(rankFrom, 2*rankFrom, getSlice((rankFrom+1)*N/size-1), N*N);
    requests[1] = world.irecv(rankFrom, 2*rankFrom, getSlice((rankFrom+1)*N/size-1), N*N);
    // printf("1\n");
  }
  // Receive first x slice of next portion
  if (rankFrom == (rank+1) % size){
    // MPI_Irecv(getSlice(rankFrom*N/size), N*N*sizeof(Vertex), MPI_BYTE, rankFrom, rankFrom, MPI_COMM_WORLD, &leftRec);
    // rightRec = world.irecv(rankFrom, rankFrom, getSlice(rankFrom*N/size), N*N);
    // world.recv(rankFrom, 2*rankFrom+1, getSlice(rankFrom*N/size), N*N);
    requests[3] = world.irecv(rankFrom, 2*rankFrom+1, getSlice(rankFrom*N/size), N*N);
    // printf("3\n");
  }
}

void VertexLattice::wait(){
  // MPI_Wait(&leftSend, MPI_STATUS_IGNORE);
  // MPI_Wait(&leftRec, MPI_STATUS_IGNORE);
  // MPI_Wait(&rightSend, MPI_STATUS_IGNORE);
  // MPI_Wait(&rightRec, MPI_STATUS_IGNORE);
  // printf("Wait 1\n");
  // leftSend.wait();
  // // printf("Wait 2\n");
  // leftRec.wait();
  // // printf("Wait 3\n");
  // rightSend.wait();
  // // printf("Wait 4\n");
  // rightRec.wait();

  // printf("Bon\n");
  
  mpi::wait_all(requests, requests+4);

  // printf("Et la ?\n");
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
