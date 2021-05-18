/*!
 *  \file file.c
 *  \brief collection of file operations
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <string.h>

FILE *fopen_safe(const char *fname, const char *mode)
{
  FILE *fp = NULL;
  fp = fopen(fname, mode);
  if (fp == NULL)
  {
    fprintf(stderr, "Could not open file `%s'\n", fname);
    exit(EXIT_FAILURE);
  }
  return fp;
}

/*!
 * Reads a file and returns the number of lines in this file.
 * Empty lines and lines starts with the comment character are neglected.
 */

int countLines(const char *fname, const char comments)
{
  FILE *fp = fopen_safe(fname, "r");
  int lines = 0;
  int32_t c;
  int32_t start = ' ';
  int flag_newline = 1;
  int flag_empty = 0;

  while (EOF != (c = fgetc(fp)))
  {
    if (flag_newline && isprint(c) && !isspace(c))
    {
      start = c;
      flag_newline = 0;
      flag_empty = 0;
    }
    if (c == '\n')
    {
      if (start != comments && !flag_empty)
      {
        ++lines;
      }
      flag_newline = 1;
      flag_empty = 1;
    }
  }
  if (c != '\n')
  {
    if (start != comments && !flag_empty)
    {
      ++lines;
    }
  }
  fclose(fp);
  return lines;
}

/*!
 * Read in formated data into arrays.
 */

void read_formated_data(const char *fname, int row, int column, const char comments, int flag, const char *fmt, ...)
{
  FILE *fp = NULL;
  fp = fopen_safe(fname, "r");
  char buf[100];

  int size = 0;
  int i = 0, j = 0, l = 0, length = 0;
  int flag_newline = 1;
  int32_t start = ' ';

  double d, **data_list_0 = NULL, *data_list_1 = NULL, ***data_list_2 = NULL;
  va_list arguments;
  va_start(arguments, fmt);

  switch (flag)
  {
  case 0: // read data into several 1d array
    data_list_0 = malloc(sizeof(double *) * row);
    for (j = 0; j < column; j++)
      data_list_0[j] = va_arg(arguments, double *);
    break;
  case 1: // read data into one 1d array
    data_list_1 = va_arg(arguments, double *);
    break;
  case 2: // read data into several 2d array
    size = va_arg(arguments, int);
    data_list_2 = malloc(sizeof(double **) * row);
    for (j = 0; j < column; j++)
      data_list_2[j] = va_arg(arguments, double **);
    break;
  default:
    fprintf(stderr, "# Error: The flag to read formated data must be 0, 1, or 2\n");
    exit(EXIT_FAILURE);
    break;
  }

  i = 0; // initialize row index
  flag_newline = 1;
  while (fgets(buf, sizeof(buf), fp) != NULL)
  {
    length = strlen(buf);
    if (!feof(fp))
    {
      /* avoid incomplete elements appearing in the buf */
      l = length - 1;
      /* find the last space character */
      while (l >= 0 && !isspace(buf[l]))
        l--;
      if (l == -1)
      {
        fprintf(stderr, "# Error: The element size of file %s is too large.\n", fname);
        exit(EXIT_FAILURE);
      }
      /* reset the file pointer */
      fseek(fp, l + 1 - length, SEEK_CUR);
      /* reset the string length */
      buf[l + 1] = '\0';
      length = strlen(buf);
    }

    if (flag_newline)
    {
      l = 0;
      /* find the first nonspace character */
      while (l < length && isspace(buf[l]))
        l++;
      if (l == length) // skip empty line
        continue;
      else // record the first nonspece character
        start = buf[l];
    }

    if (start != comments) // skip comment lines
    {
      if (flag_newline)
      {
        if (i == row)
        {
          fprintf(stderr, "# Error: The format of file %s is incorrect. "
                          "Too many rows, expected %d rows.\n",
                  fname, row);
          exit(EXIT_FAILURE);
        }
        j = 0; // reset the column number
      }

      l = 0;
      while (1)
      {
        /* Thanks Han-Qing Wang for this part. */
        while (l < length && isspace(buf[l]))
          l++;
        if (l == length)
          break;
        if (j == column)
        {
          fprintf(stderr, "# Error: The format of file %s is incorrect. "
                          "Too many columns in row %d, expected %d columns.\n",
                  fname, i + 1, column);
          exit(EXIT_FAILURE);
        }

        if (sscanf(buf + l, "%lf", &d))
        {
          switch (flag)
          {
          case 0:
            data_list_0[j][i] = d;
            break;
          case 1:
            data_list_1[i * column + j] = d;
            break;
          case 2:
            data_list_2[j][i / size][i % size] = d;
            break;
          default:
            break;
          }
          j++; // increment column number
        }

        while (l < length && !isspace(buf[l]))
        {
          l++;
        }
      }
    }

    if (buf[length - 1] == '\n' || feof(fp)) // new line next time
    {
      flag_newline = 1;
      if (start != comments)
      {
        i++; // increment row number
        if (j != column)
        {
          fprintf(stderr, "# Error: The format of file %s is incorrect. "
                          "Too few columns in row %d, expected %d columns.\n",
                  fname, i, column);
          exit(EXIT_FAILURE);
        }
      }
    }
    else // the same line next time
      flag_newline = 0;
  }

  if (i != row)
  {
    fprintf(stderr, "# Error: The format of file %s is incorrect. "
                    "Too few rows, expected %d rows.\n",
            fname, row);
    exit(EXIT_FAILURE);
  }

  fclose(fp);

  switch (flag)
  {
  case 0:
    free(data_list_0);
    break;
  case 1:
    break;
  case 2:
    free(data_list_2);
    break;
  default:
    break;
  }

  return;
}

/*!
 * Put out arrays as formated data.
 */

void write_formated_data(const char *fname, int row, int column, const char comments, const char *header, int flag, const char *fmt, ...)
{
  FILE *fp = NULL;
  fp = fopen_safe(fname, "w");

  int size = 0;
  int i = 0, j = 0, l1 = 0, l2 = 0, length = 0;
  double d, **data_list_0 = NULL, *data_list_1 = NULL, ***data_list_2 = NULL;
  va_list arguments;
  va_start(arguments, fmt);

  switch (flag)
  {
  case 0:
    data_list_0 = malloc(sizeof(double *) * row);
    for (j = 0; j < column; j++)
      data_list_0[j] = va_arg(arguments, double *);
    break;
  case 1:
    data_list_1 = va_arg(arguments, double *);
    break;
  case 2:
    size = va_arg(arguments, int);
    data_list_2 = malloc(sizeof(double **) * row);
    for (j = 0; j < column; j++)
      data_list_2[j] = va_arg(arguments, double **);
    break;
  default:
    fprintf(stderr, "# Error: The flag to read formated data must be 0, 1, or 2\n");
    exit(EXIT_FAILURE);
    break;
  }

  if (header != NULL)
  {
    length = strlen(header);
    l1 = 0, l2 = 0;
    while (l2 < length)
    {
      while (l2 < length && header[l2] != '\n')
        l2++;
      fprintf(fp, "%c ", comments);
      if (header[l2] != '\n')
      {
        fwrite(header + l1, l2 - l1, 1, fp);
        fprintf(fp, "\n");
        break;
      }
      fwrite(header + l1, l2 - l1 + 1, 1, fp);
      l2++;
      l1 = l2;
    }
  }

  for (i = 0; i < row; i++)
  {
    for (j = 0; j < column; j++)
    {
      switch (flag)
      {
      case 0:
        d = data_list_0[j][i];
        break;
      case 1:
        d = data_list_1[i * column + j];
        break;
      case 2:
        d = data_list_2[j][i / size][i % size];
        break;
      default:
        break;
      }
      fprintf(fp, fmt, d);
      fprintf(fp, "\t");
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  switch (flag)
  {
  case 0:
    free(data_list_0);
    break;
  case 1:
    break;
  case 2:
    free(data_list_2);
    break;
  default:
    break;
  }
}