/* ----------------------------------------------------------------------
   AppSurfaceRxn: simple on-lattice surface reaction KMC for SPPARKS

   This app is patterned after app_ising.{h,cpp} in SPPARKS (08 Oct 2025 build),
   i.e., it:
     - sets ninteger/ndouble and calls create_arrays() in the constructor
     - binds pointers in grow_app() via iarray[0]
     - calls solve->update() after each KMC event (like AppIsing)

   Per-site integer field:
     i1 (aka "site") stores occupancy/state:
       0 = EMPTY (*)
       1 = A* 
       2 = B*

   Pair reactions (nearest-neighbor, unordered on the bond, but counted once
   per undirected edge using an ownership rule):
     r1: A* + *  -> A* + A*        rate k1
     r2: A* + B* -> B* + B*        rate k2
     r3: B* + *  -> * + * + B(g)   rate k3  (increments a counter)

   Ownership rule:
     For each neighbor pair (i,j), only process if i < j. This prevents double
     counting while still allowing unordered chemistry on the bond.

------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(surface/rxn,AppSurfaceRxn)

#else

#ifndef SPK_APP_SURFACE_RXN_H
#define SPK_APP_SURFACE_RXN_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppSurfaceRxn : public AppLattice {
 public:
  AppSurfaceRxn(class SPPARKS *, int, char **);
  ~AppSurfaceRxn() override;

  void grow_app() override;
  void init_app() override;

  double site_energy(int) override { return 0.0; }
  void site_event_rejection(int, class RandomPark *) override {}

  double site_propensity(int) override;
  void site_event(int, class RandomPark *) override;

 private:
  // parameters
  double k1_, k2_, k3_;

  // counters
  bigint n_r1_, n_r2_, n_r3_, n_bg_;

  // app-local pointer into framework per-site arrays
  int *occ_;     // points to iarray[0] (i1/"site")

  // scratch list for solver updates (same pattern as app_ising)
  int *sites_;   // list of solver-site indices to update

  enum : int { EMPTY = 0, A_STAR = 1, B_STAR = 2 };

  inline int  state(int i) const      { return occ_[i]; }
  inline void set_state(int i, int s) { occ_[i] = s; }

  inline bool owns_edge(int i, int j) const { return i < j; }

  inline bool match_r1(int si, int sj) const { return (si==A_STAR && sj==EMPTY) || (si==EMPTY && sj==A_STAR); }
  inline bool match_r2(int si, int sj) const { return (si==A_STAR && sj==B_STAR) || (si==B_STAR && sj==A_STAR); }
  inline bool match_r3(int si, int sj) const { return (si==B_STAR && sj==EMPTY) || (si==EMPTY && sj==B_STAR); }

  void apply_r1(int i, int j, int si, int sj);
  void apply_r2(int i, int j, int si, int sj);
  void apply_r3(int i, int j, int si, int sj);

  inline void maybe_add_update(int lattice_i, int &nsites);
};

}

#endif
#endif
