#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "trains.h"

/*!
 * warp find_dict for convenience
 */
static char *find_config_param(const char *key)
{
  char *value = dict_find(key);
  if (value == NULL)
  {
    fprintf(stderr, "# Error in file %s: %s is not specicied.\n", parset.param_file_pt, key);
    exit(EXIT_FAILURE);
  }
  else
    return value;
}

/*!
 * read parameter set from parameter file.
 */
void read_parset_pt()
{
  char fname[1024];
  if (thistask == roottask)
  {
    sprintf(fname, "%s/%s", parset.file_dir, parset.param_file_pt);
    if (dict_construct(fname) == 0)
      exit(EXIT_FAILURE);

    parset_pt.flag_method = atoi(find_config_param("FlagMethod"));
    check_param_range("FlagMethod", parset_pt.flag_method, 0, 2);

    if (parset.flag_rec == 0 || parset_pt.flag_method == 0)
    {
      parset_pt.flag_evolve = atoi(find_config_param("FlagEvolve"));
      check_param_range("FlagEvolve", parset_pt.flag_evolve, 0, 1);
    }
    else
    {
      parset_pt.flag_evolve = 0;
    }

    strcpy(parset_pt.ps_catolog, find_config_param("PulsarCatalog"));
    sprintf(fname, "%s/%s", parset.file_dir, parset_pt.ps_catolog);
    parset_pt.Np = countLines(fname, '#');

    if (parset.flag_rec)
    {
      parset_pt.flag_method = atoi(find_config_param("FlagMethod"));
      check_param_range("FlagMethod", parset_pt.flag_method, 0, 2);
      strcpy(parset_pt.ptr_file, find_config_param("PTRFile"));
      strcpy(parset_pt.ptr_recon_file, find_config_param("PTRConstructFileOut"));
      if (parset_pt.flag_method == 0)
        parset_pt.Ns = atoi(find_config_param("SourceNumber"));
      else
        parset_pt.Ns = 1;
      sprintf(fname, "%s/%s", parset.file_dir, parset_pt.ptr_file);
      parset_pt.Nt = countLines(fname, '#') / parset_pt.Np;
    }
    else
    {
      parset_pt.Ns = atoi(find_config_param("SourceNumber"));
      parset_pt.Nt = atoi(find_config_param("TimingNumber"));
    }

    strcpy(parset_pt.str_par_fix, find_config_param("SourceParFix"));
    strcpy(parset_pt.str_par_fix_val, find_config_param("SourceParFixVal"));

    sprintf(fname, "%s/%s", parset.file_dir, "data/param_pt.txt");
    dict_free(fname, parset.param_file_pt);
  }
  MPI_Bcast(&parset_pt, sizeof(parset_pt), MPI_BYTE, roottask, MPI_COMM_WORLD);
  return;
}

/*!
 * read dataset.
*/
void read_data_pt()
{
  int i = 0;
  char fname[1024];

  alloc_data_pt();
  if (thistask == roottask)
  {
    sprintf(fname, "%s/%s", parset.file_dir, parset_pt.ps_catolog);
    read_formated_data(fname, parset_pt.Np, sizeof(PSRmodel) / sizeof(double), '#', 1, "%lf", psr_data);
    // write_formated_data("test.txt", parset_pt.Np, 4, '#', "test", 1, "%lf", psr_data);
    if (parset.flag_rec)
    {
      sprintf(fname, "%s/%s", parset.file_dir, parset_pt.ptr_file);
      read_formated_data(fname, parset_pt.Np * parset_pt.Nt, 2, '#', 2, "%lf", parset_pt.Nt, psr_time_data, psr_res_data);
    }
  }
  MPI_Bcast(psr_data, sizeof(PSRmodel) * parset_pt.Np, MPI_BYTE, roottask, MPI_COMM_WORLD);
  if (parset.flag_rec)
  {
    for (i = 0; i < parset_pt.Np; i++)
    {
      MPI_Bcast(psr_time_data[i], parset_pt.Nt, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
      MPI_Bcast(psr_res_data[i], parset_pt.Nt, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    }
  }

  return;
}

void alloc_data_pt()
{
  int i = 0;
  psr_data = malloc(parset_pt.Np * sizeof(PSRmodel));
  if (parset.flag_rec)
  {
    psr_time_data = malloc(parset_pt.Np * sizeof(double *));
    psr_res_data = malloc(parset_pt.Np * sizeof(double *));
    for (i = 0; i < parset_pt.Np; i++)
    {
      psr_time_data[i] = malloc(parset_pt.Nt * sizeof(double));
      psr_res_data[i] = malloc(parset_pt.Nt * sizeof(double));
    }
  }
}

void free_data_pt()
{
  int i = 0;
  free(psr_data);
  if (parset.flag_rec)
  {
    for (i = 0; i < parset_pt.Np; i++)
    {
      free(psr_time_data[i]);
      free(psr_res_data[i]);
    }
    free(psr_time_data);
    free(psr_res_data);
  }
  return;
}
