/*!
 *  \file init.c
 *  \brief initialize the program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "trains.h"

/*! 
 * This function initialize the program.
 */
void init()
{
  init_num_params();
  init_func();
  allocate_memory();
  set_par_range_model();
  set_par_fix_model();
  set_par_prior_model();

  /* initialize GSL */
  gsl_T = gsl_rng_default;
  gsl_r = gsl_rng_alloc(gsl_T);

#ifndef Debug
  if (parset.flag_rng_seed != 1)
  {
    gsl_rng_set(gsl_r, time(NULL) + thistask + 1350);
  }
  else
  {
    gsl_rng_set(gsl_r, parset.rng_seed + thistask + 1350);
  }
#else
  if (parset.flag_rng_seed != 1)
  {
    gsl_rng_set(gsl_r, 6666 + thistask + 1350);
    printf("# debugging, task %d brains random seed %d.\n", thistask, 6666 + thistask + 1350);
  }
  else
  {
    gsl_rng_set(gsl_r, parset.rng_seed + thistask + 1350);
  }
#endif
}

/*! 
 * This function initialize number of parameters.
 */
void init_num_params()
{
  int i, s, *ptr_pos, *ptr_num;

  parnum.source = parset_pt.Ns * sizeof(GWsource) / sizeof(double);
  if (parset.flag_rec && parset_pt.flag_method == 0)
  {
    parnum.phase = parset_pt.Ns * parset_pt.Np;
    if (parset_pt.flag_evolve == 1)
      parnum.dis = parset_pt.Np;
    else
      parnum.dis = 0;
  }
  else
  {
    parnum.phase = 0;
    parnum.dis = 0;
  }

  ptr_num = (int *)&parnum;
  ptr_pos = (int *)&parpos;
  s = sizeof(parpos) / sizeof(int);
  ptr_pos[0] = 0;
  num_params = ptr_num[0];
  for (i = 1; i < s; i++)
  {
    ptr_pos[i] = ptr_pos[i - 1] + ptr_num[i - 1];
    num_params += ptr_num[i];
  }
  return;
}

/*! 
 * This function initialize function handles for different BLR models.
 */
void init_func()
{
  return;
}

/*!
 * This function allocates memory for variables used throughout the code. 
 */
void allocate_memory()
{
  int i = 0;
  par_fix = malloc(num_params * sizeof(int));
  par_fix_val = malloc(num_params * sizeof(double));
  par_prior_model = malloc(num_params * sizeof(int));
  par_range_model = malloc(num_params * sizeof(double *));
  par_prior_gaussian = malloc(num_params * sizeof(double *));
  for (i = 0; i < num_params; i++)
  {
    par_range_model[i] = malloc(2 * sizeof(double));
    par_prior_gaussian[i] = malloc(2 * sizeof(double));
  }
  return;
}

/*! 
 * This function free memory.
 */
void free_memory()
{
  int i = 0;
  free(par_fix);
  free(par_fix_val);
  free(par_prior_model);
  for (i = 0; i < num_params; i++)
  {
    free(par_range_model[i]);
    free(par_prior_gaussian[i]);
  }
  free(par_range_model);
  free(par_prior_gaussian);
  gsl_rng_free(gsl_r);
  return;
}

/*!
 * This function set range of all parameters.
 */
void set_par_range_model()
{
  int i = 0, np;
  np = sizeof(GWsource) / sizeof(double);
  for (i = 0; i < parset_pt.Ns; i++)
  {
    set_source_range_model(par_range_model + parpos.source + i * np);
  }
  for (i = 0; i < parnum.phase; i++)
  {
    par_range_model[parpos.phase + i][0] = 0;
    par_range_model[parpos.phase + i][1] = M_PI;
  }
  for (i = 0; i < parnum.dis; i++)
  {
    par_range_model[parpos.dis + i][0] = 1e2;
    par_range_model[parpos.dis + i][1] = 1e4;
  }
  return;
}

/*!
 * This function copes with parameter fixing.
 */
void set_par_fix_model()
{
  int i = 0;
  int np, idx_tm, jt;

  for (i = 0; i < num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }

  set_par_fix_str(par_fix + parpos.source, par_fix_val + parpos.source,
                  parset_pt.str_par_fix, parset_pt.str_par_fix_val, parnum.source, "GW source");

  if (!parset_pt.flag_evolve)
  {
    np = sizeof(GWsource) / sizeof(double);
    idx_tm = offsetof(GWsource, log_tm) / sizeof(double);
    for (i = 0; i < parset_pt.Ns; i++)
    {
      jt = parpos.source + i * np + idx_tm;
      par_fix[jt] = 1;
      par_fix_val[jt] = log10(1e8);
    }
  }
  return;
}

/*
 *  This function set prior distribution of all parameters.
 */
void set_par_prior_model()
{
  int i = 0;
  for (i = 0; i < parpos.dis; i++)
  {
    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }
  for (i = parpos.dis; i < num_params; i++)
  {
    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = psr_data[i - parpos.dis].d;
    par_prior_gaussian[i][1] = psr_data[i - parpos.dis].delta_d;
  }
  return;
}

/*!
 * This function copes with parameter fixing using strings in param file.
 */
void set_par_fix_str(int *label, double *value, const char *label_str, const char *value_str, int num, const char *tag)
{
  int i = 0;
  int length = strlen(label_str);
  const char *pstr = value_str;

  if (thistask == roottask)
  {
    for (i = 0; i < length; i++)
    {
      if (label_str[i] == '1')
      {
        if (pstr == NULL)
        {
          printf("# %d-th %s parameter value is not provided (counting from 0).\n", i, tag);
          exit(EXIT_FAILURE);
        }
        label[i] = 1;
        sscanf(pstr, "%lf", value + i);
        printf("# %d-th %s parameter fixed, value= %f.\n", i, tag, value[i]);
        pstr = strchr(pstr, ':'); /* values are separated by ":" */
        if (pstr != NULL)
        {
          pstr++;
        }
      }
    }
  }
  MPI_Bcast(label, num, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(value, num, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  return;
}

UPDATE_FLAG get_update_flag(int pos)
{
  if (pos < 0)
    return UPDATE_ALL;
  if (pos >= parpos.source && pos < parpos.source + parnum.source)
    return UPDATE_SOURCE;
  if (pos >= parpos.phase && pos < num_params)
    return UPDATE_PULSAR;
  return UPDATE_NONE;
}
