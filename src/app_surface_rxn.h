/* ----------------------------------------------------------------------
   AppSurfaceRxn: on-lattice surface chemistry KMC for SPPARKS

   Generic surface-chemistry interface (via input commands):
     species <name> <id>           # map surface species name to integer id
     gas <name> <conc>             # define gas species concentration
     well_mixed yes|no             # shuffle lattice after each event (serial only)
     diffusion <species> <rate> [gas <name>]
                                  # A* + * -> * + A*
     adsorption <species> <rate> [gas <name>]
                                  # * -> A*
     desorption <species> <rate> [gas <name>]
                                  # A* -> *
     reaction <r1> <r2> -> <p1> <p2> <rate> [gas <name>]
                                  # adjacency-required surface reaction

   Gas dependencies:
     If "gas <name>" is provided for an event, its rate is
       k = k0 * c_gas
     where c_gas is the concentration defined by the gas command.

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

#include <map>
#include <string>
#include <vector>

namespace SPPARKS_NS {

class AppSurfaceRxn : public AppLattice {
 public:
  AppSurfaceRxn(class SPPARKS *, int, char **);
  ~AppSurfaceRxn() override;

  void grow_app() override;
  void init_app() override;
  void input_app(char *, int, char **) override;

  double site_energy(int) override { return 0.0; }
  void site_event_rejection(int, class RandomPark *) override {}

  double site_propensity(int) override;
  void site_event(int, class RandomPark *) override;

 private:
  enum RxnKind : int { RXN_UNARY = 0, RXN_PAIR = 1 };

  struct Reaction {
    int kind;     // RXN_UNARY or RXN_PAIR
    int r1, r2;   // reactants (r2 unused for unary)
    int p1, p2;   // products  (p2 unused for unary)
    double k0;    // base rate
    double rate;  // effective rate (includes gas dependency)
    int gas;      // gas index for dependency, or -1
  };

  struct Candidate {
    int j;            // neighbor index for pair reactions, -1 for unary
    int rxn;          // index into unary_rxns_ or pair_rxns_
    double r;         // event rate
    unsigned char kind; // RXN_UNARY or RXN_PAIR
    unsigned char flip; // 1 if pair matched in reverse order
  };

  // generic reaction bookkeeping
  std::vector<Reaction> unary_rxns_;
  std::vector<Reaction> pair_rxns_;
  std::vector<Candidate> candidates_;

  std::map<std::string,int> surf_index_;
  std::map<std::string,int> gas_index_;
  std::vector<std::string> gas_name_;
  std::vector<double> gas_conc_;

  int max_state_;
  bool has_user_events_;
  bool has_custom_defs_;
  int cand_max_;
  bool well_mixed_;
  std::vector<int> mix_sites_;

  // app-local pointer into framework per-site arrays
  int *occ_;     // points to iarray[0] (i1/"site")

  // scratch list for solver updates (same pattern as app_ising)
  int *sites_;   // list of solver-site indices to update

  enum : int { EMPTY = 0 };

  inline int  state(int i) const      { return occ_[i]; }
  inline void set_state(int i, int s) { occ_[i] = s; }

  inline bool owns_edge(int i, int j) const { return i < j; }

  inline int match_pair(const Reaction &rxn, int si, int sj) const {
    if (si == rxn.r1 && sj == rxn.r2) return 0;
    if (si == rxn.r2 && sj == rxn.r1) return 1;
    return -1;
  }

  inline void maybe_add_update(int lattice_i, int &nsites);

  // input helpers
  void add_species(int, char **);
  void set_gas(int, char **);
  void add_diffusion(int, char **);
  void add_adsorption(int, char **);
  void add_desorption(int, char **);
  void add_reaction(int, char **);

  int parse_surface_state(const char *);
  int parse_gas_index(const char *);
  void parse_rate_and_gas(int, char **, int, double &, int &);

  void add_unary_reaction(int r1, int p1, double k0, int gas);
  void add_pair_reaction(int r1, int r2, int p1, int p2, double k0, int gas);
  void update_effective_rates();
  void mix_after_event(class RandomPark *);
  inline bool using_custom() const { return has_user_events_; }
};

}

#endif
#endif
