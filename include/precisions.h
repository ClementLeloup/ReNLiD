// precisions.h

#include "mpi.h"

#ifndef _PRECISIONS_H_
#define _PRECISIONS_H_

// In case we want to modify float precision
typedef double FloatType;
typedef unsigned int uint;

// Spatial axes
enum Axis{x, y, z};

const int N = 64; // size of the box
extern int size; // number of processes
extern int rank;

#endif
