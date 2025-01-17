// lattice.h

#ifndef _LATTICE_H_
#define _LATTICE_H_

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
namespace mpi = boost::mpi;

#include "precisions.h"
#include "scalar.h"
#include "group.h"




class Vertex // To be modified to include multiple fields, and multi-dimensional fields
{

  friend class boost::serialization::access;
  template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & phi;
    }
  
public:
  Vertex(){}
  Vertex(RealScalar psi) : phi(psi){}//{phi=psi;}
  RealScalar getRealScalar();
  void setRealScalar(RealScalar psi);
  Vertex operator+(Vertex v);
  
protected:
  RealScalar phi;
};


// Populate edges
class Edge : public SU2
{
public:
  Edge() : SU2(){}
  Edge(const SU2& A) : SU2(A){}
};

// Indices of lattice vertices
class VertexIndex
{
public:
  VertexIndex(uint i, uint j, uint k);
  VertexIndex(const VertexIndex& vi);

  // The lattice is sliced along x
  uint getIndx() const; // get the x index
  uint getIndyz() const; // get the index in the lattice corresponding to (y,z) in a given x slice

  // Move in the periodic lattice
  VertexIndex operator+(Axis a) const;
  VertexIndex operator-(Axis a) const;
  
protected:
  uint ind[3]; // the lattice is 3-dimensional
  
};


// Indices of lattice edges
class EdgeIndex : public VertexIndex
{
public:
  EdgeIndex(uint i, uint j, uint k, Axis a) : VertexIndex(i,j,k), dir(a) {}
  EdgeIndex(const VertexIndex& vi, Axis a) : VertexIndex(vi), dir(a) {}

  uint getSlice() const;
  uint getDir() const;

protected:
  Axis dir;
};




// Vertices of the lattice
class VertexLattice
{
public:
  VertexLattice();
  VertexLattice(const VertexLattice& vl);
  VertexLattice(const Vertex* v);
  ~VertexLattice();
  Vertex& operator()(const VertexIndex& ind); // get the vertex at index ind from the lattice
  const Vertex& operator()(const VertexIndex& ind) const;
  VertexLattice& operator=(const VertexLattice& vl);
  void init();

  // MPI communication processes
  void send(int rankTo);
  void receive(int rankFrom);
  void wait();

  friend class Evolver;
  
protected:
  void alloc();
  void free();
  Vertex* getSlice(uint ind0);
  const Vertex* getSlice(uint ind0) const;

  Vertex* vertices;
  mpi::request requests[4];
  
};


// Edges of the lattice
class EdgeLattice
{
public:
  EdgeLattice();
  EdgeLattice(const EdgeLattice& el);
  ~EdgeLattice();
  Edge& operator()(const EdgeIndex& ind);
  const Edge& operator()(const EdgeIndex& ind) const;
  const Edge& operator()(const EdgeIndex& ind, Axis a) const;
  EdgeLattice& operator=(const EdgeLattice& el);
  static void goToNext(EdgeLattice& el1, EdgeLattice& el2, EdgeLattice& el3);
  void init();

  // MPI communication processes
  void send(int rank);
  void receive(int rank);
  void wait();
  
protected:
  void alloc();
  void free();
  Edge* getSlice(uint ind0, Axis a);
  const Edge* getSlice(uint ind0, Axis a) const;
  
  Edge* edges;
};



#endif
