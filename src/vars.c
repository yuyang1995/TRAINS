#include "vars.h"

/* variables for MPICH */
int thistask, totaltask, namelen;
int roottask;
char proc_name[MPI_MAX_PROCESSOR_NAME];

/* variables for GSL */
const gsl_rng_type *gsl_T;
gsl_rng *gsl_r;

/* variables for DNest */
char dnest_options_file[TRAINS_MAX_STR_LENGTH];
int which_parameter_update, which_particle_update; // which parameter and particle to be updated
int which_level_update;
double logz;

/* variables for program parameters */
PARSET parset;
PARNUM parnum;
PARPOS parpos;
int num_params;
double **par_range_model;
int *par_prior_model;
double **par_prior_gaussian;
int *par_fix, npar_fix;
double *par_fix_val;

/* Units */
double TimingUnit;
double MergeUnit;

void *best_model, *best_model_std;