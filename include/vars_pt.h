#ifndef _VARS_PT_H
#define _VARS_PT_H

#include <stdio.h>
#include "model.h"
#include "parset.h"

/* variables for program parameters */
extern PARSET_PT parset_pt;

/* data */
extern PSRmodel *psr_data;
extern double **psr_time_data, **psr_res_data;
extern double *psr_phase, **psr_time, **psr_res;
extern double *prob_pt_array, **prob_pt_particles, **prob_pt_particles_perturb;

#endif
