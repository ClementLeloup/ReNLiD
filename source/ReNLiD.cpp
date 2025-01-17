// ReNLiD.cpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <math.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <vector>
#include <fftw3.h>
#include <string>
using namespace std::complex_literals;

#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
namespace mpi = boost::mpi;

#include "progressbar.hpp"
#include "lattice.h"
#include "precisions.h"
#include "timeslice.h"
#include "evolver.h"
#include <complex>



int rank;
int size;
MPI_Comm comm;
mpi::communicator world;

// m = 0 for now
std::complex<FloatType> uk(FloatType k){
  // std::complex<FloatType> uk = (k != 0.0) ? 1.0/sqrt(2.0*k) : 0.0;
  std::complex<FloatType> uk = sqrt(lamb/(2.0*sqrt(3.0+k*k)));
  return uk;
}

// m = 0 for now
std::complex<FloatType> duk(FloatType k){
  // std::complex<FloatType> duk = (k != 0.0) ? -1i*k/sqrt(2.0*k) : 0.0;
  std::complex<FloatType> duk = -1i*sqrt(3.0+k*k + 1.0)*uk(k);
  return duk;
}




int main(int argc, char **argv){
  
  mpi::environment env;
  
  comm = world;
  rank = world.rank();
  size = world.size();

  int status = 0;  
  Evolver ev;

  std::vector<FloatType> k(N, 0.0);
  for(int i=0;i<N;i++){
    // k[i] = (2.0*M_PI*i)/N;
    k[i] = M_PI*i/(dx*N); // Because the field is real
  }

  {
    printf("Producing initial time slice\n");
    TimeSlice ts;
    fftw_complex* ak = ts.aOperator();
    std::vector<FloatType> phi(1, 1.0);
    std::vector<FloatType> dphi(1, 0.0);
    ts.init(RealScalar(phi), RealScalar(dphi), ak, &uk, &duk);
    fftw_free(ak);
    
    // printf("Before initialization energy is %f\n", ev.getEnergy());
    ev.TimeSliceToLattice(ts, 0);
  }
  
  // printf("Step %d energy = %f\n", ev.getStep(), ev.getEnergy());
  // int nstep = 200001;
  int nstep = 100;
  progressbar bar(nstep);
  bar.set_todo_char(" ");
  bar.set_done_char("█");
  // for(int s = 0; s<200001; s++){
  for(int s = 0; s<nstep; s++){  
    if(rank == 0){bar.update();}
    // save average quantities (energy, vev, ...)
    if (rank == 0 && !(s % 10)) {
      // printf("step %d\n", s);
      std::ofstream file_avg;
      file_avg.open("data/test_avg_optim_initphi.txt", std::ios::out | std::ios::app);\
      file_avg << s << "\t" << std::setprecision(15) << ev.getEnergy() << "\t" << ev.getInflaton() << "\t" << ev.getdInflaton() << "\n";
      file_avg.close();
    }

    // save power-spectrum
    if (!(s % 20000)) {
      // printf("time %f: \t Computing power-spectrum\n", s*0.005);
      TimeSlice ts_pk(ev.LatticeToTimeSlice());
      if(rank == 0){
	std::vector<FloatType> Pk = ts_pk.PowerSpectrum(k);
     
	std::ofstream file_pk;
	file_pk.open("data/test_evolver_Pk.txt", std::ios::out | std::ios::app); \
	for(int i=0;i<k.size();i++){
	  file_pk << k[i] << "\t" << Pk[i] << "\n";
	}
	file_pk.close();
      }
    }
    
    ev.evolveLattice();
    status = ev.goToNext();

    if(rank == 0 && s == (nstep - 1)){printf("\n");}
  }

  return 0;
  
}
