#ifndef MODEL_H
#define MODEL_H

typedef struct
{
  int source;
  int phase;
  int dis;
} PARNUM;

typedef struct
{
  int source;
  int phase;
  int dis;
} PARPOS;

enum UPDATE_FLAG
{
  UPDATE_ALL,
  UPDATE_SOURCE,
  UPDATE_PSR_PHASE,
  UPDATE_PSR_DIS,
  UPDATE_NONE
};
typedef enum UPDATE_FLAG UPDATE_FLAG;

typedef struct
{
  double alpha;
  double sin_delta;
  double cos_iota;
  double psi;
  double phi0;
  double log_zeta;
  double log_omega;
  double log_tm;
} GWsource;

typedef struct
{
  double alpha;
  double delta;
  double d;
  double delta_d;
  double sd;
} PSRmodel;

#endif