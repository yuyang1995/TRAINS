#include <stdio.h>

/* version number */
#define TRAINS_MAJOR_VERSION 0
#define TRAINS_MINOR_VERSION 1
#define TRAINS_PATCH_VERSION 1

/*
 * print on screen
 */
void print_version()
{
  printf("\n");
  printf("\e[1;35m"
         "%-14s: %d.%d.%d\n"
         "\e[0m",
         "TRAINS Version", TRAINS_MAJOR_VERSION, TRAINS_MINOR_VERSION, TRAINS_PATCH_VERSION);
#ifdef GITVERSION
  printf("\e[1;35m"
         "%-14s: %s\n"
         "\e[0m",
         "git log", GITVERSION);
#endif
#ifdef GITDATE
  printf("\e[1;35m"
         "%-14s: %s\n"
         "\e[0m",
         "git date", GITDATE);
#endif
  printf("\e[1;35m"
         "%-14s: %s %s\n"
         "\e[0m",
         "compiling date", __DATE__, __TIME__);
  printf("\n");
  printf("Yu-Yang Songsheng, songshengyuyang@ihep.ac.cn\n");
  return;
}

/*
 * print to file
 */
void fprint_version(FILE *fp)
{
  fprintf(fp, "# %-14s: %d.%d.%d\n", "TRAINS Version", TRAINS_MAJOR_VERSION, TRAINS_MINOR_VERSION, TRAINS_PATCH_VERSION);
#ifdef GITVERSION
  fprintf(fp, "# %-14s: %s\n", "git log", GITVERSION);
#endif
#ifdef GITDATE
  fprintf(fp, "# %-14s: %s\n", "git date", GITDATE);
#endif
  fprintf(fp, "# %-14s: %s %s\n", "compiling date", __DATE__, __TIME__);

  return;
}
