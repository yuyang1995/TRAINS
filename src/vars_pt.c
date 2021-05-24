#include "vars_pt.h"

/* variables for program parameters */
PARSET_PT parset_pt;

/* data */
PSRmodel *psr_data;
double **psr_time_data, **psr_res_data;
double *psr_phase, **psr_time, **psr_res, **psr_res_src;
double *prob_pt_array, **prob_pt_particles, **prob_pt_particles_perturb;
double ***psr_res_src_particles, ***psr_res_src_particles_perturb;
