/*!
 *  \file dnest_general.c
 *  \brief Some general function handler for dnest sampling
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dnest.h>

#include "trains.h"

/*!
 * this function generates a sample from prior.
 */
void from_prior_gen(void *model)
{
  int i;
  double *pm = (double *)model;

  for (i = 0; i < num_params; i++)
  {
    if (par_prior_model[i] == GAUSSIAN)
    {
      pm[i] = dnest_randn() * par_prior_gaussian[i][1] + par_prior_gaussian[i][0];
      dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
    }
    else
    {
      pm[i] = par_range_model[i][0] + dnest_rand() * (par_range_model[i][1] - par_range_model[i][0]);
    }
  }

  /* cope with fixed parameters */
  for (i = 0; i < num_params; i++)
  {
    if (par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }

  which_parameter_update = -1;

  return;
}

/*!
 * this function prints out parameters.
 */
void print_particle_gen(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for (i = 0; i < num_params; i++)
  {
    fprintf(fp, "%e ", pm[i]);
  }
  fprintf(fp, "\n");
  return;
}

/*!
 * This function read the particle from the file.
 */
void read_particle_gen(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for (j = 0; j < num_params; j++)
  {
    if (fscanf(fp, "%lf", psample + j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "#Error: Cannot read sample file.\n");
      exit(EXIT_FAILURE);
    }
  }
  return;
}

/*!
 * this function perturbs parameters.
 */
double perturb_gen(void *model, int which)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width;
  int which_level, size_levels;

  which_parameter_update = which;

  /* level-dependent width */
  which_level_update = dnest_get_which_level_update();
  size_levels = dnest_get_size_levels();
  which_level = which_level_update > (size_levels - 10) ? (size_levels - 10) : which_level_update;

  if (which_level > 0)
  {
    limit1 = limits[(which_level - 1) * num_params * 2 + which * 2];
    limit2 = limits[(which_level - 1) * num_params * 2 + which * 2 + 1];
    width = limit2 - limit1;
  }
  else
  {
    width = (par_range_model[which][1] - par_range_model[which][0]);
  }

  if (par_prior_model[which] == GAUSSIAN)
  {
    logH -= (-0.5 * pow((pm[which] - par_prior_gaussian[which][0]) / par_prior_gaussian[which][1], 2.0));
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5 * pow((pm[which] - par_prior_gaussian[which][0]) / par_prior_gaussian[which][1], 2.0));
  }
  else
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }

  return logH;
}

/*!
 * this function calculates likelihood.
 */
double log_likelihoods_cal_exam(const void *model)
{
  return 0.0;
}

void restart_action_gen(int iflag)
{
  return;
}

/*!
 *  print names and prior ranges for parameters 
 *
 */
void print_par_names_gen(char *fname)
{
  int i;
  FILE *fp;
  char str_fmt[TRAINS_MAX_STR_LENGTH];
  fp = fopen_safe(fname, "w");

  strcpy(str_fmt, "%4d %-30s %10.6f %10.6f %4d %4d %15.6e\n");
  printf("# Print parameter name in %s\n", fname);

  fprintf(fp, "#*************************************************\n");
  fprint_version(fp);
  fprintf(fp, "#*************************************************\n");

  fprintf(fp, "%4s %-30s %10s %10s %4s %4s %15s\n", "#", "Par", "Min", "Max", "Prior", "Fix", "Val");

  for (i = parpos.source; i < parpos.source + parnum.source; i++)
  {
    fprintf(fp, str_fmt, i, "GW Source", par_range_model[i][0], par_range_model[i][1], par_prior_model[i], par_fix[i], par_fix_val[i]);
  }

  for (i = parpos.phase; i < parpos.phase + parnum.phase; i++)
  {
    fprintf(fp, str_fmt, i, "pulsar phase", par_range_model[i][0], par_range_model[i][1], par_prior_model[i], par_fix[i], par_fix_val[i]);
  }

  for (i = parpos.dis; i < parpos.dis + parnum.dis; i++)
  {
    fprintf(fp, str_fmt, i, "pulsar distance", par_range_model[i][0], par_range_model[i][1], par_prior_model[i], par_fix[i], par_fix_val[i]);
  }
  fclose(fp);
  return;
}
