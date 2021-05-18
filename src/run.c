#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "trains.h"

/*! \file run.c
 *  \brief setup and run the program. 
 */

/*!
 * This function setups and runs the program.
 */

void begin_run()
{
  // Amplitude of timing residuals for a binary black hole with M = M_sun and DL = 1 Mpc (units: second)
  TimingUnit = pow(GRAVITY * SOLAR_MASS, 5.0 / 3.0) / CM_PER_MPC * pow(M_PI, -1.0 / 3.0) / pow(C, 4);
  // Merging time of a binary black hole with M = M_sun (units: year)
  MergeUnit = 5.0 / 256.0 * pow(C, 5) * pow(GRAVITY * SOLAR_MASS, -5.0 / 3.0) * pow(M_PI, -8.0 / 3.0) / SEC_PER_YEAR;
  /* read parameter file */
  read_parset();
  read_parset_pt();
  read_data_pt();

  init();

  if (!parset.flag_rec && thistask == roottask)
    sim();

  if (parset.flag_rec)
  {
    reconstruct_gen();
  }
  return;
}

void end_run()
{
  free_data_pt();
  free_memory();
  return;
}