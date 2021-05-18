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
    fprintf(stderr, "# Error in file %s: %s is not specicied.\n", parset.param_file, key);
    exit(EXIT_FAILURE);
  }
  else
    return value;
}

/*!
 * read parameter set from parameter file.
 */
void read_parset()
{
  char fname[1024];
  if (thistask == roottask)
  {
    sprintf(fname, "%s", parset.param_file);
    if (dict_construct(fname) == 0)
      exit(EXIT_FAILURE);

    parset.flag_rec = atoi(find_config_param("FlagRec"));
    check_param_range("FlagRec", parset.flag_rec, 0, 1);

    strcpy(parset.file_dir, find_config_param("FileDir"));
    strcpy(parset.param_file_pt, find_config_param("ParamPT"));

    sprintf(fname, "%s/%s", parset.file_dir, "data/param.txt");
    dict_free(fname, parset.param_file);
  }
  MPI_Bcast(&parset, sizeof(parset), MPI_BYTE, roottask, MPI_COMM_WORLD);
  return;
}

/*!
 * get number of particles from the option file.
 */
void get_num_particles(char *fname)
{
  FILE *fp;
  char buf[TRAINS_MAX_STR_LENGTH], buf1[TRAINS_MAX_STR_LENGTH], buf2[TRAINS_MAX_STR_LENGTH];
  fp = fopen_safe(fname, "r");

  /* default number particles */
  parset.num_particles = 1;

  while (fgets(buf, TRAINS_MAX_STR_LENGTH, fp) != NULL)
  {
    if (buf[0] == '#')
      continue;
    if (sscanf(buf, "%s%s", buf1, buf2) < 1) /* blank line */
      continue;
    if (sscanf(buf, "%s%s", buf1, buf2) < 2)
    {
      fprintf(stderr, "Error in geting number of particles.\n"
                      "Usually due to incorrect options.\n");
      exit(EXIT_FAILURE);
    }
    if (strcmp(buf1, "NumberParticles") == 0)
    {
      parset.num_particles = atoi(buf2);
      break;
    }
  }
  fclose(fp);
}

/*!
 * chech the range of flags in parameter file
 */
void check_param_range(char *name, int flag, int min, int max)
{
  if (max == INT_MAX)
  {
    if (flag < min)
    {
      fprintf(stderr, "# Error in %s: value %d is not allowed.\n"
                      "# Please specify a value greater than or equal to %d.\n",
              name, flag, min);
      exit(EXIT_FAILURE);
    }
  }
  else if (min == INT_MIN)
  {
    if (flag > max)
    {
      fprintf(stderr, "# Error in %s: value %d is not allowed.\n"
                      "# Please specify a value smaller than or equal to %d.\n",
              name, flag, max);
      exit(EXIT_FAILURE);
    }
  }
  else if (flag > max || flag < min)
  {
    fprintf(stderr, "# Error in %s: value %d is not allowed.\n"
                    "# Please specify a value in [%d, %d].\n",
            name, flag, min, max);
    exit(EXIT_FAILURE);
  }
}