// precisions.h

#include <math.h>
#include <boost/mpi.hpp>
namespace mpi = boost::mpi;

#ifndef _PRECISIONS_H_
#define _PRECISIONS_H_

// In case we want to modify float precision
typedef double FloatType;
typedef unsigned int uint;

// Spatial axes
enum Axis{x, y, z};

const int N = 64; // size of the box
const FloatType L = 16.0*M_PI;
const FloatType dtau = 0.005; // time step for the lattice evolution
const FloatType dtau_bg = 0.001; // time step for the background calculation
const FloatType dx = L/N; // lattice physical step size
const FloatType lamb = 0.0001; // scalar self-coupling

extern mpi::communicator world;
extern int size; // number of processes
extern int rank;
extern MPI_Comm comm;

// Energy and Gauss conditions' precision
const FloatType energy_precision = 0.001;
const FloatType gauss_precision = 0.001;

#endif
