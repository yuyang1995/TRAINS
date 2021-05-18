/*!
 *  \file dict.c
 *  \brief construct a dictionanry from parameter file 
 */

#include <stdio.h>
#include <stdlib.h>
#include "dict.h"

void fprint_version(FILE *fp);

static DICT *dict = NULL;

/* get value from key */
char *dict_find(const char *key)
{
  DICT *dict_item = NULL;
  HASH_FIND_STR(dict, key, dict_item);
  if (dict_item == NULL)
    return NULL;
  dict_item->value[DICT_TRAINS_MAX_STR_LENGTH - 1] = 'a'; // the item has been accessed
  return dict_item->value;
}

/* add value for key */
int dict_add(const char *key, const char *value)
{
  DICT *dict_item = NULL;
  HASH_FIND_STR(dict, key, dict_item);
  if (dict_item != NULL)
    return 0;
  else
  {
    dict_item = (DICT *)malloc(sizeof(DICT));
    strcpy(dict_item->key, key);
    strcpy(dict_item->value, value);
    dict_item->value[DICT_TRAINS_MAX_STR_LENGTH - 1] = 'n'; // the item has not been accessed
    HASH_ADD_STR(dict, key, dict_item);
    dict_item = NULL;
    return 1;
  }
}

/* construct a dictionary from a file*/
int dict_construct(const char *fname)
{
  FILE *fp;
  char str[1024], buf1[1024], buf2[1024], buf3[1024];
  dict = NULL;

  fp = fopen(fname, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    return 0;
  }
  while (fgets(str, 1024, fp) != NULL)
  {
    if (sscanf(str, "%s%s%s", buf1, buf2, buf3) < 2)
      continue;
    if (buf1[0] == '%')
      continue;

    if (dict_add(buf1, buf2) == 0)
    {
      fprintf(stderr, "# Error in file %s: Tag '%s' is multiple defined.\n", fname, buf1);
      return 0;
    }
  }
  fclose(fp);
  return 1;
}

void dict_free(const char *fname, const char *fname_par)
{
  FILE *fp;
  DICT *current_user, *tmp;

  fp = fopen(fname, "w");
  fprintf(fp, "#*************************************************\n");
  fprint_version(fp);
  fprintf(fp, "#*************************************************\n");

  fprintf(fp, "%-30s  %-s\n", "ParameterFile", fname_par);

  HASH_ITER(hh, dict, current_user, tmp)
  {
    if (current_user->value[DICT_TRAINS_MAX_STR_LENGTH - 1] == 'a')
      fprintf(fp, "%-30s  %s\n", current_user->key, current_user->value);
    HASH_DEL(dict, current_user);
    free(current_user);
  }
  fclose(fp);
}