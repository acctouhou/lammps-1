/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef MINIMIZE_CLASS
// clang-format off
MinimizeStyle(frank,MinFrank);
// clang-format on
#else

#ifndef LMP_MIN_FRANK_H
#define LMP_MIN_FRANK_H

#include "min.h"

namespace LAMMPS_NS {

class MinFrank : public Min {
 public:
  MinFrank(class LAMMPS *);
  ~MinFrank() {}
  void init();
  void setup_style();
  void reset_vectors();
  int iterate(int);

 private:
  double beta;
  
};

}    // namespace LAMMPS_NS

#endif
#endif
