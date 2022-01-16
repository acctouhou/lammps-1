// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Julien Guénolé, CNRS and
                         Erik Bitzek, FAU Erlangen-Nuernberg
------------------------------------------------------------------------- */

#include "min_frank.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "output.h"
#include "timer.h"
#include "universe.h"
#include "update.h"
#include <iostream>
#include <cmath>
#include <vector>

using namespace LAMMPS_NS;

// EPS_ENERGY = minimum normalization for energy tolerance

#define EPS_ENERGY 1.0e-8
/* ---------------------------------------------------------------------- */

MinFrank::MinFrank(LAMMPS *lmp) : Min(lmp) {}

/* ---------------------------------------------------------------------- */
double clip(double n, double lower, double upper) {
  return std::max(lower, std::min(n, upper));
}
void MinFrank::init()
{
  Min::init();

  
}

void MinFrank::setup_style()
{




}
void MinFrank::reset_vectors()
{
  nvec = 3 * atom->nlocal;
  if (nvec) xvec = atom->x[0];
  if (nvec) fvec = atom->f[0];
}

int MinFrank::iterate(int maxiter)
{
  bigint ntimestep;
  double e1=exp(1);
  double e2=exp(0.5);
  double e3=exp(-0.5);
  double e4=exp(1.05);
  double p = 0.0;
  double dfc = 0.0;
  double damp = 0.0;
  double wtf = 0.0;
  double eprevious = 0.0;
  double ecurrent = energy_force(0);
  double fdotf;
  double temp,dt,dtfm;
  double adaptive_factor;
  double **f = atom->f;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  int now_id = 0;
  int *type = atom->type;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double *m=(double *)malloc( nlocal*3 * sizeof(double) );
  double *v=(double *)malloc( nlocal*3 * sizeof(double) );
  double *s=(double *)malloc( nlocal*3 * sizeof(double) );
  double lr=0.1 ;
  double grad=0.00 ;
  double beta; 
  
  for (int j=0;j<3;j++){
		for (int i = 0; i < nlocal; i++)
		{
			now_id=i+j*nlocal;
			v[now_id]= 0.0;
			s[now_id]= 0.0;
			m[now_id]= 0.0;
		}
  }

  for (int iter = 0; iter < maxiter; iter++) {
	  
	  
	lr=0.1*cos(double(iter)/double(maxiter)*M_PI/2);
	beta=clip(1-sqrt(0.1*lr*10.0),0.0,1.0);
	double **f = atom->f;
    double **x = atom->x;
    if (timer->check_timeout(niter))
      return TIMEOUT;

    ntimestep = ++update->ntimestep;
    niter++;
	
	eprevious = ecurrent;
    ecurrent = energy_force(0);

	if (neval >= update->max_eval) return MAXEVAL;
	
	if (fabs(ecurrent-eprevious) <
        update->etol * 0.5*(fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY ))
      return ETOL;
	dt=update->dt;
	for (int j=0;j<3;j++){
		for (int i = 0; i < nlocal; i++)
		{
		now_id=i+j*nlocal;
		grad=-f[i][j]/mass[type[i]];
		p = acos(tanh(m[now_id]*grad))/M_PI;
		dfc=(1+e3)/(1+exp(-fabs(v[now_id]-p)));
		temp=grad*grad+EPS_ENERGY;
		v[now_id]=fmax(v[now_id],temp);
		damp=log(clip(e1+temp+0.5-p,e3,e4));
		adaptive_factor=lr*grad*dfc/sqrt(v[now_id]);
		m[now_id]= damp*m[now_id]*beta-adaptive_factor;
		x[i][j]+= m[now_id]-adaptive_factor;
		wtf=clip(temp/s[now_id],0.0,1.0)*fabs(p-0.5);
		v[now_id]=v[now_id]*(1.0-wtf)+wtf*temp;
		s[now_id]=temp;
		}
	}

	
	fdotf = 0.0;
    if (update->ftol > 0.0) {
      if (normstyle == MAX) fdotf = fnorm_max();        // max force norm
      else if (normstyle == INF) fdotf = fnorm_inf();   // infinite force norm
      else if (normstyle == TWO) fdotf = fnorm_sqr();   // Euclidean force 2-norm
      else error->all(FLERR,"Illegal min_modify command");
      if (fdotf < update->ftol*update->ftol) return FTOL;
    }

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }

	

  }

  return MAXITER;
}
