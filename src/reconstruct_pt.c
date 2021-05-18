/*!
 *  \file reconstruct_pt.c
 *  \brief reconstruct PTA data and recover BBH model
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_interp.h>
#include <dnest.h>

#include "trains.h"

static double prob0;

/*!
 * this function output reconstructed pulsar timing residuals
 */
void output_rec_pt(void *posterior_sample, int size_of_modeltype, int num_ps)
{
  int i, j, k;
  char fname[1024];
  void *model;
  FILE *fres;

  sprintf(fname, "%s/%s", parset.file_dir, "data/ptr_rec.txt");
  fres = fopen_safe(fname, "w");
  
  for (i = 0; i < num_ps; i++)
  {
    model = posterior_sample + i * size_of_modeltype;

    if (parset_pt.flag_method == 0)
    {
      residual_cal(model, psr_data, psr_time_data, psr_res, (double *)model + parpos.phase, (double *)model + parpos.dis, -1);
    }
    else
    {
      LLR_Mx_Av(model, psr_data, psr_time_data, psr_res_data, psr_phase);
      residual_cal(model, psr_data, psr_time_data, psr_res, psr_phase, NULL, -1);
    }

    // output pulsar timing residuals
    for (j = 0; j < parset_pt.Np; j++)
    {
      for (k = 0; k < parset_pt.Nt; k++)
      {
        fprintf(fres, "%e ", psr_res[j][k]);
      }
      fprintf(fres, "\n");
    }
  }
  fclose(fres);
}

void output_pre_pt(const void *model)
{
  char fname[1024];
  if (parset_pt.flag_method == 0)
  {
    which_particle_update = 0;
    residual_cal(model, psr_data, psr_time_data, psr_res, (double *)model + parpos.phase, (double *)model + parpos.dis, -1);
  }
  else
  {
    LLR_Mx_Av(model, psr_data, psr_time_data, psr_res_data, psr_phase);
    residual_cal(model, psr_data, psr_time_data, psr_res, psr_phase, NULL, -1);
  }

  // output timing residuals
  sprintf(fname, "%s/%s", parset.file_dir, parset_pt.ptr_recon_file);
  write_formated_data(fname, parset_pt.Np * parset_pt.Nt, 2, '#', NULL, 2, "%lf", parset_pt.Nt, psr_time_data, psr_res);
  return;
}

/*!
 * this function initializes SA reconstruction.
 */
void reconstruct_pt_init()
{
  int i;

  psr_res = malloc(parset_pt.Np * sizeof(double *));
  for (i = 0; i < parset_pt.Np; i++)
    psr_res[i] = malloc(parset_pt.Nt * sizeof(double));

  if (parset_pt.flag_method == 0)
  {
    prob_pt_particles = malloc(parset.num_particles * sizeof(double *));
    prob_pt_particles_perturb = malloc(parset.num_particles * sizeof(double *));
    for (i = 0; i < parset.num_particles; i++)
    {
      prob_pt_particles[i] = malloc(parset_pt.Np * sizeof(double));
      prob_pt_particles_perturb[i] = malloc(parset_pt.Np * sizeof(double));
    }
  }
  else
  {
    psr_phase = malloc(parset_pt.Np * sizeof(double));
    LLR_initial(parset_pt.Nt);
  }
  return;
}

/*!
 * this function finalizes SA reconstruction.
 */
void reconstruct_pt_end()
{
  int i;

  if (parset_pt.flag_method == 0)
  {
    for (i = 0; i < parset.num_particles; i++)
    {
      free(prob_pt_particles[i]);
      free(prob_pt_particles_perturb[i]);
    }
    free(prob_pt_particles);
    free(prob_pt_particles_perturb);
  }
  else
  {
    free(psr_phase);
    LLR_end();
  }

  for (i = 0; i < parset_pt.Np; i++)
  {
      free(psr_res[i]);
  }
  free(psr_res);
  return;
}

/*!
 * this function calculate probability at initial step.
 * all invoked quantities are calculated.
 */
double prob_initial_pt(const void *model)
{
  double prob = 0.0, var2, dy;
  int i, j;

  if (parset_pt.flag_method == 0)
  {
    which_particle_update = dnest_get_which_particle_update();
    which_parameter_update = -1;
    prob_pt_array = prob_pt_particles[which_particle_update];
    residual_cal(model, psr_data, psr_time_data, psr_res, (double *)model + parpos.phase, (double *)model + parpos.dis, -1);

    for (i = 0; i < parset_pt.Np; i++)
    {
      prob_pt_array[i] = 0.0;
      var2 = psr_data[i].sd * psr_data[i].sd;
      for (j = 0; j < parset_pt.Nt; j++)
      {
        dy = psr_res[i][j] - psr_res_data[i][j];
        prob_pt_array[i] += (-0.5 * (dy * dy) / var2) - 0.5 * log(var2 * 2.0 * M_PI);
      }
      prob += prob_pt_array[i];
    }
  }
  else
  {
    prob0 = 0;
    for (i = 0; i < parset_pt.Np; i++)
    {
      var2 = psr_data[i].sd * psr_data[i].sd;
      for (j = 0; j < parset_pt.Nt; j++)
      {
        dy = psr_res_data[i][j];
        prob0 += (-0.5 * (dy * dy) / var2) - 0.5 * log(var2 * 2.0 * M_PI);
      }
    }
    if (parset_pt.flag_method == 1)
      prob = LLR_Mx_Av(model, psr_data, psr_time_data, psr_res_data, psr_phase);
    else if (parset_pt.flag_method == 2)
      prob = LLR_Mx_Av(model, psr_data, psr_time_data, psr_res_data, NULL);
    prob = prob + prob0;
  }
  return prob;
}

/*!
 * this function calculates probabilities.
 *
 * At each MCMC step, only one parameter is updated, which only changes some values; thus,
 * optimization that reuses the unchanged values can improve computation efficiency.
 */
double prob_pt(const void *model)
{
  double prob = 0.0, var2, dy;
  int i, j, jp;
  UPDATE_FLAG uflag;

  if (parset_pt.flag_method == 0)
  {
    uflag = get_update_flag(which_parameter_update);
    which_particle_update = dnest_get_which_particle_update();

    if (uflag == UPDATE_SOURCE || uflag == UPDATE_ALL)
      jp = -1;
    else
      jp = (which_parameter_update - parpos.phase) % parset_pt.Np;
    residual_cal(model, psr_data, psr_time_data, psr_res, (double *)model + parpos.phase, (double *)model + parpos.dis, jp);

    prob_pt_array = prob_pt_particles_perturb[which_particle_update];
    for (i = 0; i < parset_pt.Np; i++)
    {
      if (jp == -1 || i == jp)
      {
        prob_pt_array[i] = 0.0;
	      var2 = psr_data[i].sd * psr_data[i].sd;
        for (j = 0; j < parset_pt.Nt; j++)
        {
          dy = psr_res[i][j] - psr_res_data[i][j];
          prob_pt_array[i] += (-0.5 * (dy * dy) / var2) - 0.5 * log(var2 * 2.0 * M_PI);
        }
        prob += prob_pt_array[i];
      }
      else
      {
        prob += prob_pt_particles[which_particle_update][i];
      }
    }
  }
  else if (parset_pt.flag_method == 1)
  {
    prob = LLR_Mx_Av(model, psr_data, psr_time_data, psr_res_data, psr_phase);
    prob = prob + prob0;
  }
  else if (parset_pt.flag_method == 2)
  {
    prob = LLR_Mx_Av(model, psr_data, psr_time_data, psr_res_data, NULL);
    prob = prob + prob0;
  }
  return prob;
}
