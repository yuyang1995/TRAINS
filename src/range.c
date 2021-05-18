/*!
 *  \file range.c
 *  \brief set the range of parameters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void set_source_range_model(double **range)
{
  int i = 0;
  // alpha
  range[i][0] = 0.0;
  range[i++][1] = 2 * M_PI;
  // sin(delta)
  range[i][0] = -1.0;
  range[i++][1] = 1.0;
  // cos(iota)
  range[i][0] = -1.0;
  range[i++][1] = 1.0;
  // psi
  range[i][0] = 0.0;
  range[i++][1] = M_PI;
  // phi0
  range[i][0] = 0.0;
  range[i++][1] = M_PI;
  // log(zeta)
  range[i][0] = -15.0;
  range[i++][1] = -6.0;
  // log(omega)
  range[i][0] = log10(1);
  range[i++][1] = log10(100);
  // log(tmerge)
  range[i][0] = log10(1);
  range[i++][1] = log10(1e8);
  return;
}
