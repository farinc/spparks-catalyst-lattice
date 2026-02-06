/* ----------------------------------------------------------------------
   AppSurfaceRxn implementation
------------------------------------------------------------------------- */

#include "app_surface_rxn.h"

#include "solve.h"
#include "random_park.h"
#include "error.h"

#include <cstdlib>

using namespace SPPARKS_NS;

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

void AppSurfaceRxn::init_app()
{
  delete [] sites_;
  sites_ = new int[2 + 2*maxneigh];

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (occ_[i] < 0 || occ_[i] > 2) flag = 1;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values (expected 0,1,2)");
}

double AppSurfaceRxn::site_propensity(int i)
{
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

/* ---------------------------------------------------------------------- */

void AppSurfaceRxn::apply_r1(int i, int j, int si, int sj)
{
  if (si==A_STAR && sj==EMPTY) set_state(j, A_STAR);
  else if (si==EMPTY && sj==A_STAR) set_state(i, A_STAR);
  n_r1_++;
}

void AppSurfaceRxn::apply_r2(int i, int j, int si, int sj)
{
  if (si==A_STAR && sj==B_STAR) set_state(i, B_STAR);
  else if (si==B_STAR && sj==A_STAR) set_state(j, B_STAR);
  n_r2_++;
}

void AppSurfaceRxn::apply_r3(int i, int j, int si, int sj)
{
  if (si==B_STAR && sj==EMPTY) set_state(i, EMPTY);
  else if (si==EMPTY && sj==B_STAR) set_state(j, EMPTY);
  n_r3_++;
  n_bg_++;
}

/* ---------------------------------------------------------------------- */

inline void AppSurfaceRxn::maybe_add_update(int lattice_i, int &nsites)
{
  const int isite = i2site[lattice_i];
  if (isite < 0) return;
  sites_[nsites++] = isite;
  propensity[isite] = site_propensity(lattice_i);
}

/* ----------------------------------------------------------------------
   KMC event: choose and apply one event owned by i, then update propensities
------------------------------------------------------------------------- */

void AppSurfaceRxn::site_event(int i, RandomPark *random)
{
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

  // Update propensities for endpoints and their neighbors (delpropensity=1)
  int nsites = 0;
  maybe_add_update(i, nsites);
  maybe_add_update(j, nsites);

  for (int m = 0; m < numneigh[i]; m++) maybe_add_update(neighbor[i][m], nsites);
  for (int m = 0; m < numneigh[j]; m++) maybe_add_update(neighbor[j][m], nsites);

  solve->update(nsites, sites_, propensity);
}
