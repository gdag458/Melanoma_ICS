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
AppStyle(rxn/diff,AppRxnDiff)

#else

#ifndef SPK_APP_RXN_DIFF_H
#define SPK_APP_RXN_DIFF_H

#include "app_lattice.h"

namespace SPPARKS_NS {

#define MAX_SPECIES 10

class AppRxnDiff : public AppLattice {
//  friend class DiagRxnDiff;

 public:
  AppRxnDiff(class SPPARKS *, int, char **);
  ~AppRxnDiff();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);
  

  int tcell_total_population; // total T cell population count

 protected:
  int engstyle;
  int firsttime;

  int *esites;
  int *echeck;


  int nspecies;                  // # of unique species
  char **sname;                  // ID of each species

  int nreactions;                // # of user defined reactions
  char **rname;                  // ID of each reaction
  int **localReactants;          // local reactants;
                                 // localReactants[I][J] = local number of reactants for reaction I of species J
  int **nbrReactants;            // neighbor reactants
  int **localDeltaPop;           // local population change if even is carried out;
                                 // localDeltaPop[I][J] = change in local population of species J if reaction I is carried out
  int **nbrDeltaPop;             // neighbor population change if even is carried out
  double *rate;                  // rate[I] = input rate for reaction I
  int *rxnStyle;                 // rxnStyle[I] = LOCAL or NBR for reaction I

  int *rxn_count;               // statistics

  int max_cell_area_per_voxel; // area units in a pixel
  int cytotoxic_T_area_count;  // area units of an aCD8
  int exhausted_T_area_count;  // area units of an eCD8
  int cancer_area_count;       // area units of a cancer cell
  int mac_area_count;          // area units of a TAM

  struct Event {           // one event for an owned site
    int style;             // reaction style = LOCAL, NBR
    int which;             // which reaction
    int jpartner;          // which J neighbor of I are part of event
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void clear_events(int);

  void add_event(int, int, int, double, int);
  void add_rxn(int, char **);
  void add_species(int, char **);

  int find_reaction(char *);
  int find_species(char *);

  virtual double custom_multiplier(int, int, int, int);

  void grow_reactions(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

FIXME

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.


*/
