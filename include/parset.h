#ifndef _PARSET_H
#define _PARSET_H

#define TRAINS_MAX_STR_LENGTH 256

#include "uthash.h"

typedef struct
{
  int flag_postprc;
  int flag_restart;
  int flag_sample_info;
  int flag_temp;
  double temperature;
  int flag_exam_prior;
  int flag_rng_seed, rng_seed;
  int flag_force_update;
  int flag_force_run;
  int flag_help, flag_end;
  int flag_para_name;
  int num_particles;

  int flag_rec;
  char file_dir[TRAINS_MAX_STR_LENGTH];
  char param_file[TRAINS_MAX_STR_LENGTH];
  char param_file_pt[TRAINS_MAX_STR_LENGTH];
} PARSET;

typedef struct
{
  char ps_catolog[TRAINS_MAX_STR_LENGTH];
  char ptr_file[TRAINS_MAX_STR_LENGTH];
  char ptr_recon_file[TRAINS_MAX_STR_LENGTH];
  int flag_method;
  int flag_evolve;
  int Ns, Np, Nt;
  char str_par_fix[TRAINS_MAX_STR_LENGTH];
  char str_par_fix_val[TRAINS_MAX_STR_LENGTH];
} PARSET_PT;

#endif