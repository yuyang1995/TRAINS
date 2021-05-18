#ifndef _DICT_H
#define _DICT_H

#include "uthash.h"

#define DICT_TRAINS_MAX_STR_LENGTH 256

typedef struct
{
  char key[DICT_TRAINS_MAX_STR_LENGTH]; /* key */
  char value[DICT_TRAINS_MAX_STR_LENGTH];
  UT_hash_handle hh; /* makes this structure hashable */
} DICT;

char *dict_find(const char *key);
int dict_add(const char *key, const char *value);
int dict_construct(const char *fname);
void dict_free(const char *fname, const char *fname_par);

#endif