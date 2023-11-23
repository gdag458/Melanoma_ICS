/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(rxn/diff/custom,AppRxnDiffCustom)

#else

#ifndef SPK_APP_RXN_DIFF_CUSTOM_H
#define SPK_APP_RXN_DIFF_CUSTOM_H

#include "app_rxn_diff.h"

namespace SPPARKS_NS {

class AppRxnDiffCustom : public AppRxnDiff {
 public:
  AppRxnDiffCustom(class SPPARKS *, int, char **);
  ~AppRxnDiffCustom() {}

 protected:
  int max_cell_area_per_voxel; // area in a voxel/pixel/chamber
  int cytotoxic_T_area_count;  // area taken by activated CD8+ T cells
  int exhausted_T_area_count;  // area taken by exhausted CD8+ T cells
  int cancer_area_count;       // area taken by melanoma cells
  int mac_area_count;          // area taken by TAMs
  double tcell_stick_fac;      // factor by which T cell diffusion rate changes for T cells in contact with cancer cells
  double pd1_fac;              // (not used in paper) factor which modulates rates given pd1 therapy
  int prolif_rec_thresh;       // the number of cancer cells which must be killed in last 12 hours to maintain half of max prolif/recruitment rates
  int do_carrying_capacity;    // whether or not to apply carrying capacity restrictions for CD8+ T cells
  int tcell_carrying_capacity; // carrying capacity limit for all CD8+ T cells

  double custom_multiplier(int, int, int, int);

};

}

#endif
#endif
