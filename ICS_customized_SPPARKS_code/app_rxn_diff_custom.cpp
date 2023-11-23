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

#include <iostream>

#include "math.h"
#include "string.h"
#include "app_rxn_diff_custom.h"
#include "random_park.h"
#include "error.h"
#include <random>
#include <iostream>

using namespace SPPARKS_NS;

enum{LOCAL,NBR};

/* ---------------------------------------------------------------------- */
AppRxnDiffCustom::AppRxnDiffCustom(SPPARKS *spk, int narg, char **arg) : 
  AppRxnDiff(spk,narg,arg)
{
  // parse arguments for RxnDiffCustom class only, not children
  if (strcmp(style,"rxn/diff/custom") != 0) return;
  if (narg != 11) error->all(FLERR,"Illegal app_style command"); // input arguments in in. file
  max_cell_area_per_voxel = atoi(arg[1]); // maximum area units in cell
  cancer_area_count = atoi(arg[2]);       // area units of cancer cell
  cytotoxic_T_area_count = atoi(arg[3]);  // area units of aCD8s
  exhausted_T_area_count = atof(arg[4]);  // area units of eCD8s
  mac_area_count = atof(arg[5]);          // area units of TAMs
  tcell_stick_fac = atof(arg[6]);         // factor by which diffusion is reduced for CD8s when in contact with cancer cell
  pd1_fac = atof(arg[7]);                 // factor by which PD1 blocker effects interaction rates
  prolif_rec_thresh = atof(arg[8]);       // michaelis-menten population of CD8s for half rate
  do_carrying_capacity = atof(arg[9]);    // include carrying capacity or not
  tcell_carrying_capacity = atof(arg[10]);// CD8 T cell carrying capacity
}

/* ----------------------------------------------------------------------
 * Custom rate multiplier function
------------------------------------------------------------------------- */
double AppRxnDiffCustom::custom_multiplier(int i, int rstyle, int which, int jpartner)
{
	int latCol = 100;
	int latRow = 100;


int max_filling = 4;
//data IDs to access the populations within the population array
int macID = 6;
int cd8aID = 0;
int canID = 1;
int cd8eID = 2;
int recID = 4;
int latEID = 7;
int deadID = 3;

//area units of different cell types
int maxA = max_cell_area_per_voxel;
int tcA = cytotoxic_T_area_count;
int teA = exhausted_T_area_count;
int macA = mac_area_count;
int cancA = cancer_area_count;
	
if(which == 0){ // cancer proliferation and motility
	if(which == 0){
		bool done = false;
        int neighLim = 3; // number of layers of neighbors that cancer can jump to
        int widthLim = 3 + 2*(neighLim-1);
		int width = 1;
		int ix = i%latCol;
		int iy = i/latCol;
        while(!done){
        int halfw = int((width-1) / 2);
        for(int ly=0; ly<2; ly++){ // check to see if there are any spaces to which a cancer cell may proliferate
            for(int lx=0; lx<width; lx++){
                int target1x = ix + lx - halfw;
                int target1y = iy + (2*ly - 1)*halfw;
                int target2x = ix + (2*ly-1)*halfw;
                int target2y = iy + lx - halfw;
                target1x -= floor(double(target1x) / latCol) * latCol;
                target1y -= floor(double(target1y) / latCol) * latCol;
                target2x -= floor(double(target2x)/latCol)*latCol;
                target2y -= floor(double(target2y)/latCol)*latCol;
                int target1 = target1x + latCol * target1y;
                int target2 = target2x + latCol*target2y;
                if ((population[canID][target1]*cancA + population[macID][target1]*macA < maxA)){
                    return 1.0;
                    done = true;
                }
                if(population[canID][target2]*cancA + population[macID][target2]*macA < maxA){
                    return 1.0;
                    done = true;
                }
            }
        }
        if(width >= widthLim){
            return 0.0;
            done =  true;
        }else{width += 2;}
        }
    }else{return 0.0;}

  }else if(which == 1){ // CD8 recruitment
    int Trec_dCan_denom = prolif_rec_thresh;
    int t_carry_cap = tcell_carrying_capacity;
    int t_tot_pop = tcell_total_population;
    double deadTrec_fac = (1-(double(t_tot_pop)*do_carrying_capacity/double(t_carry_cap)))*double(2*nDead/double(Trec_dCan_denom+nDead)); // Michaelis-Menten
    
    return deadTrec_fac;
  }else if ((which == 2 || which == 4 || which == 5)){//CD8 proliferation and motility
    int filling = tcA*population[cd8aID][i]+teA*population[cd8eID][i] + cancA*population[canID][i] + macA*population[macID][i];
    
    int t_carry_cap = tcell_carrying_capacity;
    int t_tot_pop = tcell_total_population; 
    int Tpro_dCan_denom = prolif_rec_thresh;

    int yOn;
    if(which == 2){ // T cell proliferation
    	if (filling + tcA > maxA) {
    		return 0.0;
    	}else{
            double deadTprolif_fac = (1-(double(t_tot_pop)*do_carrying_capacity/double(t_carry_cap)))*double(2*nDead/double(Tpro_dCan_denom+nDead))/pd1_fac; // Michaelis-Menten

		    return deadTprolif_fac;
    	}
    }else if(which == 4 || which == 5){ // Tcell motility
	  //int stickF = 0.03;
      bool problem = false;
	  if((which == 4 && population[cd8aID][i] <= 0) || (which == 5 && population[cd8eID][i] <= 0)){
	    problem = true;
	  }
      int weight;
      if(which == 4){
        weight = tcA;
      }else if(which == 5){
        weight = teA;
      }
      int nfilling = tcA*population[cd8aID][jpartner] + teA*population[cd8eID][jpartner] + cancA*population[canID][jpartner] + macA*population[macID][jpartner];
      
	  if(maxA < nfilling + weight){
	    return 0.0;
	  }else{ // heck for neighboring cancer cells to adhere
		int ix = i%latCol;
		int iy = i/latCol;
        for(int k=0;k < 2; k++){
          int targx = ix;
          int targy = iy - 1 + 2*k; 
          targy -= floor(double(targy)/latCol)*latCol;
          int targ = targx + targy*latCol;
          if(population[1][targ] > 0){
            return tcell_stick_fac;
          }
        }
        for(int h=0; h < 2; h++){
          int targx = ix - 1 + 2*h;
          int targy = iy;
          targx -= floor(double(targx)/latCol)*latCol;
          int targ = targx + targy*latCol;
          if(population[canID][targ] > 0){
            return tcell_stick_fac;
          }
        }
		return 1.0;
	  }
    }
  }else if(which == 8){ // no longer using
	return 0.0;
  }else if(which == 10){ // Macrophage motility
    int filling = tcA*population[cd8aID][i] + teA*population[cd8eID][i] + cancA*(population[canID][i] + macA*population[macID][i]);
    int nfilling = tcA*population[cd8aID][jpartner] + teA*population[cd8eID][jpartner] + cancA*(population[canID][jpartner] + macA*population[macID][jpartner]);
    if(maxA - nfilling < macA){ // block motion if neighboring chamber is full
      return 0.0;
    }else{return 1.0;}
  }
  if(which == 7 || which == 11 || which == 13 || which == 15){ // apply pd1 ratio to these rates
    return pd1_fac;
  }
  if(which == 3 || which == 14){
    return 1/pd1_fac;
  }
  
  return 1.0;
}

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

#include <iostream>

#include "math.h"
#include "string.h"
#include "app_rxn_diff_custom.h"
