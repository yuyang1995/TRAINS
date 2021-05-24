/*!
 *  \file sim.c
 *  \brief generate mock data
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "trains.h"

void *model;
double *rho;

void sim()
{
  char fname[1024];
  int i, j, k, np;
  double sigma;

  sim_init();

  double *pm = (double *)model;
  // read_formated_data("data/sim_source.txt", parset_pt.Ns, sizeof(GWsource) / sizeof(double), '#', 1, "%lf", pm);
  np = sizeof(GWsource) / sizeof(double);
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_source.txt");
  write_formated_data(fname, parset_pt.Ns, np, '#', NULL, 1, "%lf", pm);

  residual_cal(pm, psr_data, psr_time, psr_res_src, psr_phase, NULL, -1, -1);
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_phase.txt");
  write_formated_data(fname, parset_pt.Ns, parset_pt.Np, '#', NULL, 1, "%lf", psr_phase);

  for (i = 0; i < parset_pt.Ns; i++)
  {
    residual_shift(psr_res_src + i * parset_pt.Np, -1);
    rho[i] = 0;
    for (j = 0; j < parset_pt.Np; j++)
    {
      sigma = psr_data[j].sd;
      for (k = 0; k < parset_pt.Nt; k++)
        rho[i] += pow(psr_res_src[i * parset_pt.Np + j][k] / sigma, 2);
    }
    rho[i] = sqrt(rho[i]);
  }
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_snr.txt");
  write_formated_data(fname, parset_pt.Ns, 1, '#', NULL, 0, "%e", rho);

  if (parset_pt.Ns > 1)
  {
    residual_sum(psr_res_src, psr_res, -1);
  }
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_signal.txt");
  write_formated_data(fname, parset_pt.Np * parset_pt.Nt, 2, '#', NULL, 2, "%e", parset_pt.Nt, psr_time, psr_res);

  for (i = 0; i < parset_pt.Np; i++)
    for (j = 0; j < parset_pt.Nt; j++)
      psr_res[i][j] += gsl_ran_ugaussian(gsl_r) * psr_data[i].sd;
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_ptr.txt");
  write_formated_data(fname, parset_pt.Np * parset_pt.Nt, 2, '#', NULL, 2, "%e", parset_pt.Nt, psr_time, psr_res);

  sim_end();
  return;
}

void sim_init()
{
  int i, j;

  model = malloc(num_params * sizeof(double));
  double *pm = (double *)model;
  set_source_value_sim(pm + parpos.source);

  psr_phase = malloc(parset_pt.Ns * parset_pt.Np * sizeof(double));
  psr_time = malloc(parset_pt.Np * sizeof(double *));
  psr_res = malloc(parset_pt.Np * sizeof(double *));
  for (i = 0; i < parset_pt.Np; i++)
  {
    psr_time[i] = malloc(parset_pt.Nt * sizeof(double));
    psr_res[i] = malloc(parset_pt.Nt * sizeof(double));
    for (j = 0; j < parset_pt.Nt; j++)
    {
      psr_time[i][j] = 0 + j * 14 / 365.25;
    }
  }

  rho = malloc(parset_pt.Ns * sizeof(double));
  if (parset_pt.Ns == 1)
  {
    psr_res_src = psr_res;
  }
  else
  {
    psr_res_src = malloc(parset_pt.Ns * parset_pt.Np * sizeof(double *));
    for (i = 0; i < parset_pt.Ns * parset_pt.Np; i++)
    {
      psr_res_src[i] = malloc(parset_pt.Nt * sizeof(double));
    }
  }

  for (i = 0; i < num_params; i++)
  {
    if (par_fix[i] == 1)
    {
      pm[i] = par_fix_val[i];
    }
  }

  return;
}

void sim_end()
{
  int i;
  free(model);

  free(psr_phase);
  for (i = 0; i < parset_pt.Np; i++)
  {
    free(psr_time[i]);
    free(psr_res[i]);
  }
  free(psr_time);
  free(psr_res);

  if (parset_pt.Ns > 1)
  {
    for (i = 0; i < parset_pt.Ns * parset_pt.Np; i++)
    {
      free(psr_res_src[i]);
    }
    free(psr_res_src);
  }
  free(rho);

  return;
}

/* 
 * set parameter values for simulation
 */
void set_source_value_sim(double *pm)
{
  int i, jz, jo, jt;
  int np = sizeof(GWsource) / sizeof(double);
  int idx_zeta = offsetof(GWsource, log_zeta) / sizeof(double);
  int idx_omega = offsetof(GWsource, log_omega) / sizeof(double);
  int idx_tm = offsetof(GWsource, log_tm) / sizeof(double);
  double fgw, ds, Mc, Amp, tmerge;

  for (i = 0; i < num_params; i++)
  {
    pm[i] = par_range_model[i][0] + gsl_rng_uniform(gsl_r) * (par_range_model[i][1] - par_range_model[i][0]);
  }

  for (i = 0; i < parset_pt.Ns; i++)
  {
    jz = parpos.source + i * np + idx_zeta;
    jo = parpos.source + i * np + idx_omega;
    jt = parpos.source + i * np + idx_tm;
    if (parset_pt.flag_evolve == 0)
    {
      // from log_omega to fgw
      fgw = pow(10, pm[jo]) / 2 / M_PI / SEC_PER_YEAR;
      // distance of GW source, [100, 1000] Mpc
      ds = pow(1e-3 + (1 - 1e-3) * gsl_rng_uniform(gsl_r), 1.0 / 3.0) * 1000;
      // generate log uniform distribution of mass, [1e6, 1e10] Msun
      Mc = pow(10, 6 + 4 * gsl_rng_uniform(gsl_r));
      // the amplitude of timing residuals
      Amp = TimingUnit / ds * pow(Mc, 5.0 / 3.0) * pow(fgw, -1.0 / 3.0);
      pm[jz] = log10(Amp);
      tmerge = MergeUnit * pow(Mc, -5.0 / 3.0) * pow(fgw, -8.0 / 3.0);
      pm[jt] = log10(tmerge);
    }
    else
    {
      pm[jz] = -15 + 8 * gsl_rng_uniform(gsl_r);
      pm[jt] = 2 + 6 * gsl_rng_uniform(gsl_r);
    }
  }

  return;
}
