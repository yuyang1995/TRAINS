#ifndef _VARS_H
#define _VARS_H

#include <float.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_const_cgsm.h>

#include "model.h"
#include "parset.h"

#define GRAVITY GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT
#define SOLAR_MASS GSL_CONST_CGSM_SOLAR_MASS
#define C GSL_CONST_CGSM_SPEED_OF_LIGHT
#define SEC_PER_YEAR 31557600.0
#define CM_PER_PC GSL_CONST_CGSM_PARSEC
#define CM_PER_MPC (1.0e6 * CM_PER_PC)

#define EPS (DBL_MIN) /* epsilon of the machine as in Matlab */

enum PRIOR_TYPE
{
    GAUSSIAN = 1,
    UNIFORM = 2
};

/* variables for MPICH */
extern int thistask, totaltask, namelen;
extern int roottask;
extern char proc_name[MPI_MAX_PROCESSOR_NAME];

/* variables for GSL */
extern const gsl_rng_type *gsl_T;
extern gsl_rng *gsl_r;

/* variables for DNest */
extern char dnest_options_file[TRAINS_MAX_STR_LENGTH];
extern int which_parameter_update, which_particle_update; // which parameter and particle to be updated
extern int which_level_update;
extern double *limits; // external from dnest
extern double logz;

/* variables for program parameters */
extern PARSET parset;
extern PARNUM parnum;
extern PARPOS parpos;
extern int num_params;
extern double **par_range_model;
extern int *par_prior_model;
extern double **par_prior_gaussian;
extern int *par_fix, npar_fix;
extern double *par_fix_val;
extern int num_params_source;
extern double **source_range_model;

/* Units */
extern double TimingUnit;
extern double MergeUnit;

extern void *best_model, *best_model_std;
#endif