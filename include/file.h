#ifndef _FILE_H
#define _FILE_H

#include <stdio.h>

FILE *fopen_safe(const char *fname, const char *mode);
int countLines(const char *fname, const char comment);
void read_formated_data(const char *fname, int row, int column, const char comments, int flag, const char *fmt, ...);
void write_formated_data(const char *fname, int row, int column, const char comments, const char *header, int flag, const char *fmt, ...);
#endif