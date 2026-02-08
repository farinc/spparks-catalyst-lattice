/* ----------------------------------------------------------------------
   AppSurfaceRxn implementation
------------------------------------------------------------------------- */

#include "app_surface_rxn.h"

#include "solve.h"
#include "random_park.h"
#include "error.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace SPPARKS_NS;

namespace {

bool is_number(const char *str)
{
  if (str == nullptr || *str == '\0') return false;
  char *end = nullptr;
  std::strtod(str,&end);
  return end != str && *end == '\0';
}

bool is_empty_name(const std::string &name)
{
  if (name == "*") return true;
  if (name.size() != 5) return false;
  return (std::toupper(static_cast<unsigned char>(name[0])) == 'E' &&
          std::toupper(static_cast<unsigned char>(name[1])) == 'M' &&
          std::toupper(static_cast<unsigned char>(name[2])) == 'P' &&
          std::toupper(static_cast<unsigned char>(name[3])) == 'T' &&
          std::toupper(static_cast<unsigned char>(name[4])) == 'Y');
}

}

AppSurfaceRxn::AppSurfaceRxn(SPPARKS *spk, int narg, char **arg)
  : AppLattice(spk,narg,arg)
{
  ninteger = 1;
  ndouble  = 0;

  delpropensity = 1;
  delevent      = 1;

  allow_kmc       = 1;
  allow_rejection = 0;
  allow_masking   = 0;

  create_arrays();

  k1_ = k2_ = k3_ = 1.0;
  if (narg == 4) {
    k1_ = atof(arg[1]);
    k2_ = atof(arg[2]);
    k3_ = atof(arg[3]);
  } else if (narg != 1) error->all(FLERR,"Illegal app_style surface/rxn [k1 k2 k3]");

  max_state_ = EMPTY;
  has_user_events_ = false;
  has_custom_defs_ = false;
  cand_max_ = 0;
  well_mixed_ = false;

  surf_index_["EMPTY"] = EMPTY;

  occ_ = NULL;
  sites_ = NULL;
}

void AppSurfaceRxn::grow_app()
{
  occ_ = iarray[0];
}

AppSurfaceRxn::~AppSurfaceRxn()
{
  delete [] sites_;
}

void AppSurfaceRxn::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"species") == 0) add_species(narg,arg);
  else if (strcmp(command,"gas") == 0) set_gas(narg,arg);
  else if (strcmp(command,"diffusion") == 0 || strcmp(command,"diffuse") == 0)
    add_diffusion(narg,arg);
  else if (strcmp(command,"adsorption") == 0) add_adsorption(narg,arg);
  else if (strcmp(command,"desorption") == 0) add_desorption(narg,arg);
  else if (strcmp(command,"reaction") == 0) add_reaction(narg,arg);
  else if (strcmp(command,"well_mixed") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal well_mixed command");
    if (strcmp(arg[0],"yes") == 0) well_mixed_ = true;
    else if (strcmp(arg[0],"no") == 0) well_mixed_ = false;
    else error->all(FLERR,"Illegal well_mixed command");
  }
  else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   input helpers
------------------------------------------------------------------------- */

void AppSurfaceRxn::add_species(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal species command");
  if (!is_number(arg[1])) error->all(FLERR,"Illegal species command");

  std::string name(arg[0]);
  if (!name.empty() && name.back() == '*') name.pop_back();
  if (name.empty()) error->all(FLERR,"Illegal species command");

  const int id = atoi(arg[1]);
  if (id < 0) error->all(FLERR,"Illegal species command");

  if (is_empty_name(name) && id != EMPTY)
    error->all(FLERR,"EMPTY must map to state 0");

  auto it = surf_index_.find(name);
  if (it != surf_index_.end()) {
    if (it->second != id) error->all(FLERR,"Species ID already defined");
    return;
  }

  for (auto &kv : surf_index_) {
    if (kv.second == id) error->all(FLERR,"Species ID already defined");
  }

  surf_index_[name] = id;
  max_state_ = std::max(max_state_, id);
  has_custom_defs_ = true;
}

void AppSurfaceRxn::set_gas(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal gas command");
  if (!is_number(arg[1])) error->all(FLERR,"Illegal gas command");

  const double conc = atof(arg[1]);
  if (conc < 0.0) error->all(FLERR,"Illegal gas command");

  std::string name(arg[0]);
  if (name.empty()) error->all(FLERR,"Illegal gas command");

  auto it = gas_index_.find(name);
  if (it == gas_index_.end()) {
    const int idx = static_cast<int>(gas_conc_.size());
    gas_index_[name] = idx;
    gas_name_.push_back(name);
    gas_conc_.push_back(conc);
  } else {
    gas_conc_[it->second] = conc;
  }

  has_custom_defs_ = true;
}

int AppSurfaceRxn::parse_surface_state(const char *token)
{
  if (token == nullptr || *token == '\0')
    error->all(FLERR,"Unknown surface species in surface/rxn command");

  if (strcmp(token,"*") == 0) return EMPTY;

  std::string name(token);
  if (!name.empty() && name.back() == '*') name.pop_back();
  if (is_empty_name(name)) return EMPTY;

  if (is_number(name.c_str())) {
    const int id = atoi(name.c_str());
    if (id < 0) error->all(FLERR,"Unknown surface species in surface/rxn command");
    return id;
  }

  auto it = surf_index_.find(name);
  if (it == surf_index_.end())
    error->all(FLERR,"Unknown surface species in surface/rxn command");
  return it->second;
}

int AppSurfaceRxn::parse_gas_index(const char *token)
{
  if (token == nullptr || *token == '\0')
    error->all(FLERR,"Unknown gas species in surface/rxn command");
  auto it = gas_index_.find(token);
  if (it == gas_index_.end())
    error->all(FLERR,"Unknown gas species in surface/rxn command");
  return it->second;
}

void AppSurfaceRxn::parse_rate_and_gas(int narg, char **arg, int start,
                                       double &k0, int &gas)
{
  int i = start;
  if (i >= narg) error->all(FLERR,"Reaction has no rate");

  if (strcmp(arg[i],"rate") == 0) i++;
  if (i >= narg) error->all(FLERR,"Reaction has no rate");
  if (!is_number(arg[i])) error->all(FLERR,"Reaction has no numeric rate");

  k0 = atof(arg[i]);
  if (k0 < 0.0) error->all(FLERR,"Reaction has no numeric rate");
  i++;

  gas = -1;
  if (i < narg) {
    if (strcmp(arg[i],"gas") != 0) error->all(FLERR,"Illegal reaction command");
    if (i+1 >= narg) error->all(FLERR,"Illegal reaction command");
    gas = parse_gas_index(arg[i+1]);
    i += 2;
  }

  if (i != narg) error->all(FLERR,"Illegal reaction command");
}

void AppSurfaceRxn::add_unary_reaction(int r1, int p1, double k0, int gas)
{
  Reaction rxn;
  rxn.kind = RXN_UNARY;
  rxn.r1 = r1;
  rxn.r2 = -1;
  rxn.p1 = p1;
  rxn.p2 = -1;
  rxn.k0 = k0;
  rxn.rate = k0;
  rxn.gas = gas;
  unary_rxns_.push_back(rxn);
  has_user_events_ = true;
  max_state_ = std::max(max_state_, std::max(r1,p1));
}

void AppSurfaceRxn::add_pair_reaction(int r1, int r2, int p1, int p2,
                                      double k0, int gas)
{
  Reaction rxn;
  rxn.kind = RXN_PAIR;
  rxn.r1 = r1;
  rxn.r2 = r2;
  rxn.p1 = p1;
  rxn.p2 = p2;
  rxn.k0 = k0;
  rxn.rate = k0;
  rxn.gas = gas;
  pair_rxns_.push_back(rxn);
  has_user_events_ = true;
  max_state_ = std::max(max_state_, std::max(std::max(r1,r2),std::max(p1,p2)));
}

void AppSurfaceRxn::add_diffusion(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal diffusion command");
  const int species = parse_surface_state(arg[0]);
  double k0;
  int gas;
  parse_rate_and_gas(narg,arg,1,k0,gas);
  add_pair_reaction(species,EMPTY,EMPTY,species,k0,gas);
}

void AppSurfaceRxn::add_adsorption(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal adsorption command");
  const int product = parse_surface_state(arg[0]);
  double k0;
  int gas;
  parse_rate_and_gas(narg,arg,1,k0,gas);
  add_unary_reaction(EMPTY,product,k0,gas);
}

void AppSurfaceRxn::add_desorption(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal desorption command");
  const int reactant = parse_surface_state(arg[0]);
  double k0;
  int gas;
  parse_rate_and_gas(narg,arg,1,k0,gas);
  add_unary_reaction(reactant,EMPTY,k0,gas);
}

void AppSurfaceRxn::add_reaction(int narg, char **arg)
{
  if (narg < 6) error->all(FLERR,"Illegal reaction command");

  const int r1 = parse_surface_state(arg[0]);
  const int r2 = parse_surface_state(arg[1]);
  if (strcmp(arg[2],"->") != 0) error->all(FLERR,"Illegal reaction command");
  const int p1 = parse_surface_state(arg[3]);
  const int p2 = parse_surface_state(arg[4]);

  double k0;
  int gas;
  parse_rate_and_gas(narg,arg,5,k0,gas);
  add_pair_reaction(r1,r2,p1,p2,k0,gas);
}

void AppSurfaceRxn::update_effective_rates()
{
  for (size_t i = 0; i < unary_rxns_.size(); i++) {
    Reaction &rxn = unary_rxns_[i];
    rxn.rate = rxn.k0;
    if (rxn.gas >= 0) {
      if (rxn.gas >= (int)gas_conc_.size())
        error->all(FLERR,"Unknown gas species in surface/rxn command");
      rxn.rate *= gas_conc_[rxn.gas];
    }
  }

  for (size_t i = 0; i < pair_rxns_.size(); i++) {
    Reaction &rxn = pair_rxns_[i];
    rxn.rate = rxn.k0;
    if (rxn.gas >= 0) {
      if (rxn.gas >= (int)gas_conc_.size())
        error->all(FLERR,"Unknown gas species in surface/rxn command");
      rxn.rate *= gas_conc_[rxn.gas];
    }
  }

  cand_max_ = static_cast<int>(unary_rxns_.size()) +
              static_cast<int>(pair_rxns_.size()) * maxneigh;
  if (cand_max_ < 1) cand_max_ = 1;
  candidates_.resize(cand_max_);
}

void AppSurfaceRxn::mix_after_event(RandomPark *random)
{
  if (!well_mixed_) return;
  if (nset != 1) error->all(FLERR,"well_mixed requires no sectors");
  if (nprocs > 1) error->all(FLERR,"well_mixed requires a single processor");

  for (int i = nlocal - 1; i > 0; --i) {
    const int j = static_cast<int>(random->uniform() * (i + 1));
    const int tmp = occ_[i];
    occ_[i] = occ_[j];
    occ_[j] = tmp;
  }

  const int nsites = set[0].nlocal;
  if ((int)mix_sites_.size() != nsites) {
    mix_sites_.resize(nsites);
    for (int i = 0; i < nsites; i++) mix_sites_[i] = i;
  }

  int *site2i = set[0].site2i;
  for (int isite = 0; isite < nsites; isite++) {
    const int lattice_i = site2i[isite];
    propensity[isite] = site_propensity(lattice_i);
  }

  solve->update(nsites, mix_sites_.data(), propensity);
}

void AppSurfaceRxn::init_app()
{
  delete [] sites_;
  sites_ = new int[2 + 2*maxneigh];

  if (!using_custom()) {
    if (has_custom_defs_)
      error->all(FLERR,"No reactions defined for surface/rxn app");

    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if (occ_[i] < 0 || occ_[i] > 2) flag = 1;

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall) error->all(FLERR,"One or more sites have invalid values (expected 0,1,2)");
    return;
  }

  update_effective_rates();

  if (well_mixed_) {
    if (nset != 1) error->all(FLERR,"well_mixed requires no sectors");
    if (nprocs > 1) error->all(FLERR,"well_mixed requires a single processor");
    mix_sites_.resize(nlocal);
    for (int i = 0; i < nlocal; i++) mix_sites_[i] = i;
  }

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (occ_[i] < 0 || occ_[i] > max_state_) flag = 1;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) {
    char msg[128];
    sprintf(msg,"One or more sites have invalid values (expected 0..%d)",max_state_);
    error->all(FLERR,msg);
  }
}

double AppSurfaceRxn::site_propensity(int i)
{
  if (!using_custom()) {
    double a = 0.0;
    const int si = state(i);

    for (int m = 0; m < numneigh[i]; m++) {
      const int j = neighbor[i][m];
      if (!owns_edge(i,j)) continue;

      const int sj = state(j);
      if (match_r1(si,sj)) a += k1_;
      if (match_r2(si,sj)) a += k2_;
      if (match_r3(si,sj)) a += k3_;
    }
    return a;
  }

  double a = 0.0;
  const int si = state(i);

  for (size_t r = 0; r < unary_rxns_.size(); r++) {
    const Reaction &rxn = unary_rxns_[r];
    if (si == rxn.r1) a += rxn.rate;
  }

  for (int m = 0; m < numneigh[i]; m++) {
    const int j = neighbor[i][m];
    if (!owns_edge(i,j)) continue;

    const int sj = state(j);
    for (size_t r = 0; r < pair_rxns_.size(); r++) {
      const Reaction &rxn = pair_rxns_[r];
      if (match_pair(rxn,si,sj) >= 0) a += rxn.rate;
    }
  }

  return a;
}

/* ---------------------------------------------------------------------- */

void AppSurfaceRxn::apply_r1(int i, int j, int si, int sj)
{
  if (si==A_STAR && sj==EMPTY) set_state(j, A_STAR);
  else if (si==EMPTY && sj==A_STAR) set_state(i, A_STAR);
}

void AppSurfaceRxn::apply_r2(int i, int j, int si, int sj)
{
  if (si==A_STAR && sj==B_STAR) set_state(i, B_STAR);
  else if (si==B_STAR && sj==A_STAR) set_state(j, B_STAR);
}

void AppSurfaceRxn::apply_r3(int i, int j, int si, int sj)
{
  if (si==B_STAR && sj==EMPTY) set_state(i, EMPTY);
  else if (si==EMPTY && sj==B_STAR) set_state(j, EMPTY);
}

/* ---------------------------------------------------------------------- */

inline void AppSurfaceRxn::maybe_add_update(int lattice_i, int &nsites)
{
  const int isite = i2site[lattice_i];
  if (isite < 0) return;
  sites_[nsites++] = isite;
  propensity[isite] = site_propensity(lattice_i);
}

void AppSurfaceRxn::site_event(int i, RandomPark *random)
{
  if (!using_custom()) {
    struct Cand { int j; int rxn; double r; int si; int sj; };
    Cand cand[256];
    int nc = 0;
    double sum = 0.0;

    const int si = state(i);

    for (int m = 0; m < numneigh[i]; m++) {
      const int j = neighbor[i][m];
      if (!owns_edge(i,j)) continue;

      const int sj = state(j);

      if (match_r1(si,sj)) { cand[nc++] = {j, 1, k1_, si, sj}; sum += k1_; }
      if (match_r2(si,sj)) { cand[nc++] = {j, 2, k2_, si, sj}; sum += k2_; }
      if (match_r3(si,sj)) { cand[nc++] = {j, 3, k3_, si, sj}; sum += k3_; }

      if (nc >= 256) error->all(FLERR,"AppSurfaceRxn: too many candidates; increase buffer");
    }

    if (sum <= 0.0) return;

    double u = random->uniform() * sum;
    double c = 0.0;
    int pick = nc-1;
    for (int k = 0; k < nc; k++) {
      c += cand[k].r;
      if (u < c) { pick = k; break; }
    }

    const int j   = cand[pick].j;
    const int rxn = cand[pick].rxn;
    const int si0 = cand[pick].si;
    const int sj0 = cand[pick].sj;

    if (rxn == 1) apply_r1(i, j, si0, sj0);
    else if (rxn == 2) apply_r2(i, j, si0, sj0);
    else if (rxn == 3) apply_r3(i, j, si0, sj0);

    if (well_mixed_) {
      mix_after_event(random);
      return;
    }

    // Update propensities for endpoints and their neighbors (delpropensity=1)
    int nsites = 0;
    maybe_add_update(i, nsites);
    maybe_add_update(j, nsites);

    for (int m = 0; m < numneigh[i]; m++) maybe_add_update(neighbor[i][m], nsites);
    for (int m = 0; m < numneigh[j]; m++) maybe_add_update(neighbor[j][m], nsites);

    solve->update(nsites, sites_, propensity);
    return;
  }

  int nc = 0;
  double sum = 0.0;
  const int si = state(i);

  for (size_t r = 0; r < unary_rxns_.size(); r++) {
    const Reaction &rxn = unary_rxns_[r];
    if (si == rxn.r1) {
      if (nc >= cand_max_)
        error->all(FLERR,"AppSurfaceRxn: too many candidates; increase buffer");
      candidates_[nc++] = { -1, static_cast<int>(r), rxn.rate,
                            static_cast<unsigned char>(RXN_UNARY), 0 };
      sum += rxn.rate;
    }
  }

  for (int m = 0; m < numneigh[i]; m++) {
    const int j = neighbor[i][m];
    if (!owns_edge(i,j)) continue;

    const int sj = state(j);
    for (size_t r = 0; r < pair_rxns_.size(); r++) {
      const Reaction &rxn = pair_rxns_[r];
      const int flip = match_pair(rxn,si,sj);
      if (flip >= 0) {
        if (nc >= cand_max_)
          error->all(FLERR,"AppSurfaceRxn: too many candidates; increase buffer");
        candidates_[nc++] = { j, static_cast<int>(r), rxn.rate,
                              static_cast<unsigned char>(RXN_PAIR),
                              static_cast<unsigned char>(flip) };
        sum += rxn.rate;
      }
    }
  }

  if (sum <= 0.0) return;

  double u = random->uniform() * sum;
  double c = 0.0;
  int pick = nc-1;
  for (int k = 0; k < nc; k++) {
    c += candidates_[k].r;
    if (u < c) { pick = k; break; }
  }

  const Candidate &cand = candidates_[pick];

  if (cand.kind == RXN_UNARY) {
    const Reaction &rxn = unary_rxns_[cand.rxn];
    set_state(i, rxn.p1);
  } else {
    const Reaction &rxn = pair_rxns_[cand.rxn];
    const int j = cand.j;
    if (cand.flip == 0) {
      set_state(i, rxn.p1);
      set_state(j, rxn.p2);
    } else {
      set_state(i, rxn.p2);
      set_state(j, rxn.p1);
    }
  }

  if (well_mixed_) {
    mix_after_event(random);
    return;
  }

  // Update propensities for endpoints and their neighbors (delpropensity=1)
  int nsites = 0;
  maybe_add_update(i, nsites);
  if (cand.kind == RXN_PAIR) maybe_add_update(cand.j, nsites);

  for (int m = 0; m < numneigh[i]; m++) maybe_add_update(neighbor[i][m], nsites);
  if (cand.kind == RXN_PAIR) {
    for (int m = 0; m < numneigh[cand.j]; m++)
      maybe_add_update(neighbor[cand.j][m], nsites);
  }

  solve->update(nsites, sites_, propensity);
}
