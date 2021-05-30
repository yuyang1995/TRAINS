/*!
 *  \file range.c
 *  \brief set the range of parameters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "trains.h"

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

void output_source_range_model(double *mean, double *std)
{
  int i = 0, j = 0;
  FILE *fp;
  char range_file[TRAINS_MAX_STR_LENGTH];
  double ll, ul;
  sprintf(range_file, "%s/%s", parset.file_dir, "data/range_pt.txt");
  fp = fopen_safe(range_file, "w");
  for (i = parpos.source; i < parpos.source + parnum.source; i++)
  {
    j = (i - parpos.source) % num_params_source;
    ll = source_range_model[j][0];
    ul = source_range_model[j][1];
    if (par_fix[i] == 0)
    {
      if (ll < mean[i] - 5 * std[i])
        ll = mean[i] - 5 * std[i];
      if (ul > mean[i] + 5 * std[i])
        ul = mean[i] + 5 * std[i];
    }
    fprintf(fp, "%f %f\n", ll, ul);
  }
  fclose(fp);
  return;
}

void input_source_range_model(double **range, int size)
{
  int i = 0;
  char range_file[TRAINS_MAX_STR_LENGTH];
  sprintf(range_file, "%s/%s", parset.file_dir, "data/range_pt.txt");
  double *data_buf;
  data_buf = malloc(size * 2 * sizeof(double));

  read_formated_data(range_file, size, 2, '#', 1, "%lf", data_buf);
  for (i = 0; i < size; i++)
  {
    range[i][0] = data_buf[i * 2 + 0];
    range[i][1] = data_buf[i * 2 + 1];
  }
  free(data_buf);
  return;
}

double log_prior_compress_rate()
{
  int i, j;
  double r, logR = 0;
  for (i = parpos.source; i < parpos.source + parnum.source; i++)
  {
    j = (i - parpos.source) % num_params_source;
    r = (par_range_model[i][1] - par_range_model[i][0]) / (source_range_model[j][1] - source_range_model[j][0]);
    logR += log(r);
  }
  return logR;
}