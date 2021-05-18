#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dnest.h>

#include "trains.h"

/*!
 *  calculate mean and standard deviation of a sample.
 */
void sample_stats(double *sample, double *mean, double *std, int num_ps)
{
  int i, j;

  for (j = 0; j < num_params; j++)
    mean[j] = std[j] = 0.0;

  for (i = 0; i < num_ps; i++)
    for (j = 0; j < num_params; j++)
      mean[j] += sample[i * num_params + j];

  for (j = 0; j < num_params; j++)
    mean[j] /= num_ps;

  for (i = 0; i < num_ps; i++)
    for (j = 0; j < num_params; j++)
      std[j] += pow(sample[i * num_params + j] - mean[j], 2.0);

  for (j = 0; j < num_params; j++)
  {
    if (num_ps > 1)
      std[j] = sqrt(std[j] / (num_ps - 1.0));
    else
      std[j] = 0.0;
  }
  for (j = 0; j < num_params; j++)
    printf("Best params %d %f +- %f\n", j, mean[j], std[j]);
  return;
}

/*!
 * postprocessing.
 */
void postprocess_gen()
{
  char posterior_sample_file[TRAINS_MAX_STR_LENGTH];
  int num_ps;
  void *posterior_sample;
  int size_of_modeltype = num_params * sizeof(double);
  FILE *fp;

  best_model = malloc(size_of_modeltype);
  best_model_std = malloc(size_of_modeltype);

  if (thistask == roottask)
  {
    // get number of lines in posterior sample file
    dnest_get_posterior_sample_file(posterior_sample_file);
    fp = fopen_safe(posterior_sample_file, "r");
    if (fscanf(fp, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(EXIT_FAILURE);
    }
    printf("# Number of points in posterior sample: %d\n", num_ps);
    fclose(fp);
    // read posterior sample
    posterior_sample = malloc(num_ps * size_of_modeltype);
    read_formated_data(posterior_sample_file, num_ps, num_params, '#', 1, "%lf", posterior_sample);
    
    output_rec_pt(posterior_sample, size_of_modeltype, num_ps);

    sample_stats(posterior_sample, best_model, best_model_std, num_ps);
    free(posterior_sample);
  }
  return;
}

/*!
 * this function run dnest sampleing, 
 * reconstruct line profile and differential phase curve using the best estimates for parameters.
 */
void reconstruct_gen()
{
  int i;

  char restart_fname[1024], tag[10];

  sprintf(tag, "_pt");
  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "options/OPTIONSPT");
  sprintf(restart_fname, "/data/restart%s_dnest.txt", tag);

  //configure restart of dnest
  int argc = 0;
  char **argv;
  argv = malloc(11 * sizeof(char *));
  for (i = 0; i < 11; i++)
  {
    argv[i] = malloc(TRAINS_MAX_STR_LENGTH * sizeof(char));
  }
  //setup argc and argv
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc], parset.file_dir);
  strcat(argv[argc++], restart_fname);

  if (parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], restart_fname);
  }

  if (parset.flag_postprc == 1)
  {
    strcpy(argv[argc++], "-p");
  }
  if (parset.flag_temp == 1)
  {
    sprintf(argv[argc++], "-t%f", parset.temperature);
  }
  if (parset.flag_sample_info == 1)
  {
    strcpy(argv[argc++], "-c");
  }

  //level-dependent sampling
  strcpy(argv[argc++], "-l");

  // sample tag
  strcpy(argv[argc++], "-g");
  strcpy(argv[argc++], tag);

  // get number of particles
  if (thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  reconstruct_pt_init();
  dnest_pt(argc, argv);
  

  if (parset.flag_exam_prior != 1 && parset.flag_para_name != 1)
  {
    postprocess_gen();
    if (thistask == roottask)
    {
      output_pre_pt(best_model);
    }
  }

  reconstruct_pt_end();
  
  for (i = 0; i < 11; i++)
  {
    free(argv[i]);
  }
  free(argv);

  if (thistask == roottask)
  {
    printf("Ends reconstruct%s.\n", tag);
  }
  return;
}