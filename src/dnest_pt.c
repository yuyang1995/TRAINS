/*!
 *  \file dnest_pt.c
 *  \brief run dnest sampling for SA data analysis.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <dnest.h>

#include "trains.h"

DNestFptrSet *fptrset_pt;

double log_likelihoods_cal_initial_pt(const void *model);
double log_likelihoods_cal_restart_pt(const void *model);
double log_likelihoods_cal_pt(const void *model);
void from_prior_pt(void *model);
double perturb_pt(void *model);
void accept_action_pt();
void kill_action_pt(int i, int i_copy);

/*!
 * this function does dnest samling.
 */
int dnest_pt(int argc, char **argv)
{
  fptrset_pt = dnest_malloc_fptrset();

  fptrset_pt->from_prior = from_prior_pt;
  fptrset_pt->read_particle = read_particle_gen;
  fptrset_pt->print_particle = print_particle_gen;
  fptrset_pt->restart_action = restart_action_gen;

  fptrset_pt->accept_action = accept_action_pt;
  fptrset_pt->kill_action = kill_action_pt;
  fptrset_pt->perturb = perturb_pt;

  if (parset.flag_exam_prior != 1)
  {
    fptrset_pt->log_likelihoods_cal_initial = log_likelihoods_cal_initial_pt;
    fptrset_pt->log_likelihoods_cal_restart = log_likelihoods_cal_restart_pt;
    fptrset_pt->log_likelihoods_cal = log_likelihoods_cal_pt;
  }
  else
  {
    fptrset_pt->log_likelihoods_cal_initial = log_likelihoods_cal_exam;
    fptrset_pt->log_likelihoods_cal_restart = log_likelihoods_cal_exam;
    fptrset_pt->log_likelihoods_cal = log_likelihoods_cal_exam;
  }

  char fname[TRAINS_MAX_STR_LENGTH];
  if (thistask == roottask)
  {
    sprintf(fname, "%s/%s", parset.file_dir, "data/para_names_pt.txt");
    print_par_names_gen(fname);
  }

  if (parset.flag_para_name != 1)
    logz = dnest(argc, argv, fptrset_pt, num_params, NULL, NULL, NULL, "data/", dnest_options_file, NULL, NULL);

  dnest_free_fptrset(fptrset_pt);

  return 0;
}

/*!
 * this function generates a sample from prior.
 */
void from_prior_pt(void *model)
{
  from_prior_gen(model);
  int count = parset_pt.Ns;
  int size = sizeof(GWsource) / sizeof(double);
  int index = offsetof(GWsource, log_omega) / sizeof(double);
  gsl_sort((double *)model + index, size, count);
  return;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_initial_pt(const void *model)
{
  double logL;
  logL = prob_initial_pt(model);
  return logL;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_restart_pt(const void *model)
{
  double logL;
  logL = prob_initial_pt(model);
  return logL;
}

/*!
 * this function calculates likelihood.
 */
double log_likelihoods_cal_pt(const void *model)
{
  double logL;
  logL = prob_pt(model);
  return logL;
}

/*!
 * this function perturbs parameters.
 */
double perturb_pt(void *model)
{
  int which;
  do
  {
    which = dnest_rand_int(num_params);
  } while (par_fix[which] == 1);
  double logH = perturb_gen(model, which);

  int count = parset_pt.Ns;
  int size = sizeof(GWsource) / sizeof(double);
  int index = offsetof(GWsource, log_omega) / sizeof(double);
  gsl_sort((double *)model + index, size, count);
  return logH;
}

void accept_action_pt()
{
  int jp;
  UPDATE_FLAG uflag;

  if (parset_pt.flag_method == 0)
  {
    uflag = get_update_flag(which_parameter_update);
    if (uflag == UPDATE_SOURCE || uflag == UPDATE_ALL)
    {
      memcpy(prob_pt_particles[which_particle_update], prob_pt_particles_perturb[which_particle_update], parset_pt.Np * sizeof(double));
    }
    else
    {
      jp = (which_parameter_update - parpos.phase) % parset_pt.Np;
      prob_pt_particles[which_particle_update][jp] = prob_pt_particles_perturb[which_particle_update][jp];
    }
  }
  return;
}

/*
 * action when particle i is killed in cdnest sampling.
 * particle i_copy's properties is copyed to particle i. 
 */
void kill_action_pt(int i, int i_copy)
{
  if (parset_pt.flag_method == 0)
  {
    memcpy(prob_pt_particles[i], prob_pt_particles[i_copy], parset_pt.Np * sizeof(double));
  }
  return;
}
