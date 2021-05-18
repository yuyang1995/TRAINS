/*!
 *  \file MxAvPhase.c
 *  \brief calculate log-likelihood of PTA data using max phase and average phase method
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_double.h>

#include "trains.h"

static double *X, *Y, *Z;
static gsl_poly_complex_workspace *w1;
static gsl_integration_workspace *w2;

/*!
 * This function calculates inner product of two time series with error bars
 */
static double InnProduct(int N, double *x, double *y, double sd)
{
  int i;
  double result = 0;

  for (i = 0; i < N; i++)
    result += (x[i] * y[i]);
  result /= (sd * sd);
  return result;
}

/*!
 * log-likelihood ratio, Equation (15) in Wang et al. (2015)
 */
static double logLambda(double x, double *b)
{
  double value = b[0] + b[1] * cos(x) + b[2] * sin(x) + b[3] * sin(2.0 * x) +
                 b[4] * cos(x) * cos(x);
  return value;
}

/*!
 * likelihood ratio
 */
static double Lambda(double x, void *p)
{
  double *b = (double *)p;
  double value = exp(logLambda(x, b) - b[5]);
  return value;
}

/*!
 * This function calculates log-likelihood ratio by
 * (1) Maximization over pulsar phases
 * (2) Marginalization over pulsar phases
 */
double LLR_Mx_Av(const void *pm, PSRmodel *p, double **t, double **res_data, double *phase)
{
  int j, k;
  GWsource *s = (GWsource *)pm;
  double alphaS = s->alpha;
  double sin_deltaS = s->sin_delta;
  double cos_inc = s->cos_iota;
  double psi = s->psi;
  double phi0 = s->phi0;
  double zeta = pow(10, s->log_zeta) / 2.0;
  double omega = pow(10, s->log_omega);
  double cos_deltaS = sqrt(1 - sin_deltaS * sin_deltaS);

  double sin_deltaP, cos_deltaP, alphatilde, cos_theta;
  double Pp, Pc, Fp, Fc, A1, A2, A, phiA, omegat;
  double ip[9];
  double c[5], e[5], z[8], llh_tmp[10], phase_tmp[10];
  double b[6], result, error;
  double LLR = 0.0;

  size_t index;
  gsl_function F;
  F.function = &Lambda;
  F.params = &(b[0]);

  for (j = 0; j < parset_pt.Np; j++)
  {
    alphatilde = alphaS - p[j].alpha;
    cos_deltaP = cos(p[j].delta);
    sin_deltaP = sin(p[j].delta);

    Pp = -pow(cos_deltaP, 2) * (1 - 2 * pow(cos(alphatilde), 2) + pow(cos(alphatilde), 2) * pow(cos_deltaS, 2)) + pow(sin_deltaP, 2) * pow(cos_deltaS, 2) - 2 * (sin_deltaP * cos_deltaP) * cos(alphatilde) * (sin_deltaS * cos_deltaS);

    Pc = -2 * cos_deltaP * sin(alphatilde) * (cos_deltaP * cos(alphatilde) * sin_deltaS - sin_deltaP * cos_deltaS);

    cos_theta = cos_deltaP * cos(p[j].alpha) * cos_deltaS * cos(alphaS) + cos_deltaP * sin(p[j].alpha) * cos_deltaS * sin(alphaS) + sin_deltaP * sin_deltaS;

    Fp = Pp / (1 - cos_theta);
    Fc = Pc / (1 - cos_theta);

    A1 = (1 + cos_inc * cos_inc) * (Fp * cos(2 * psi) + Fc * sin(2 * psi));
    A2 = 2 * cos_inc * (Fp * sin(2 * psi) - Fc * cos(2 * psi));
    A = 2 * zeta * sqrt(A1 * A1 + A2 * A2);
    phiA = atan2(A1, A2);

    for (k = 0; k < parset_pt.Nt; k++)
    {
      omegat = omega * t[j][k];
      X[k] = 0.5 * A * cos(phiA + omegat);
      Y[k] = -0.5 * A * sin(phiA + omegat);
      Z[k] = -0.5 * A * cos(2 * phi0 + phiA + omegat);
    }

    ip[0] = InnProduct(parset_pt.Nt, res_data[j], X, p[j].sd);
    ip[1] = InnProduct(parset_pt.Nt, res_data[j], Y, p[j].sd);
    ip[2] = InnProduct(parset_pt.Nt, res_data[j], Z, p[j].sd);
    ip[3] = InnProduct(parset_pt.Nt, X, X, p[j].sd);
    ip[4] = InnProduct(parset_pt.Nt, X, Y, p[j].sd);
    ip[5] = InnProduct(parset_pt.Nt, X, Z, p[j].sd);
    ip[6] = InnProduct(parset_pt.Nt, Y, Y, p[j].sd);
    ip[7] = InnProduct(parset_pt.Nt, Y, Z, p[j].sd);
    ip[8] = InnProduct(parset_pt.Nt, Z, Z, p[j].sd);

    c[0] = -ip[0] + ip[5];
    c[1] = ip[1] - ip[7];
    c[2] = 0.5 * (ip[3] - ip[6]);
    c[3] = -ip[4];

    e[4] = 4 * (c[2] * c[2] + c[3] * c[3]);
    e[3] = 4 * (c[0] * c[2] + c[1] * c[3]);
    e[2] = c[0] * c[0] + c[1] * c[1] - e[4];
    e[1] = -2.0 * c[1] * c[3] - 4.0 * c[0] * c[2];
    e[0] = c[3] * c[3] - c[0] * c[0];

    gsl_poly_complex_solve(e, 5, w1, z);

    b[0] = ip[2] - 0.5 * (ip[6] + ip[8]);
    b[1] = ip[0] - ip[5];
    b[2] = ip[1] - ip[7];
    b[3] = -0.5 * ip[4];
    b[4] = 0.5 * (ip[6] - ip[3]);

    for (k = 0; k < 4; k++)
    {
      if (z[2 * k + 1] == 0.0 && fabs(z[2 * k]) <= 1.0)
      {
        phase_tmp[2 * k] = acos(z[2 * k]);
        llh_tmp[2 * k] = logLambda(phase_tmp[2 * k], b);
        phase_tmp[2 * k + 1] = -phase_tmp[2 * k];
        llh_tmp[2 * k + 1] = logLambda(phase_tmp[2 * k + 1], b);
      }
      else
      {
        phase_tmp[2 * k] = GSL_NAN;
        llh_tmp[2 * k] = GSL_NEGINF;
        phase_tmp[2 * k] = GSL_NAN;
        llh_tmp[2 * k + 1] = GSL_NEGINF;
      }
    }
    phase_tmp[8] = 0;
    llh_tmp[8] = logLambda(0, b);
    phase_tmp[9] = M_PI;
    llh_tmp[9] = logLambda(M_PI, b);

    if (phase != NULL) // maximum phase method
    {
      gsl_sort_largest_index(&index, 1, llh_tmp, 1, 10);
      phase[j] = phase_tmp[index] / 2.0;
      LLR += llh_tmp[index];
    }
    else // avarage phase method
    {
      gsl_sort_largest(b + 5, 1, llh_tmp, 1, 10);
      // gsl_integration_qags(&F, 0, 2 * M_PI, 0, 1e-5, 1000, w2, &result, &error);
      gsl_integration_qag(&F, 0, 2 * M_PI, 0, 1e-5, 1000, GSL_INTEG_GAUSS21, w2, &result, &error);
      LLR = LLR + log(result) + b[5] - log(2 * M_PI);
    }
  }
  return LLR;
}

void LLR_initial(int Nt)
{
  X = malloc(Nt * sizeof(double));
  Y = malloc(Nt * sizeof(double));
  Z = malloc(Nt * sizeof(double));
  w1 = gsl_poly_complex_workspace_alloc(5);
  w2 = gsl_integration_workspace_alloc(1000);
}

void LLR_end()
{
  free(X);
  free(Y);
  free(Z);
  gsl_poly_complex_workspace_free(w1);
  gsl_integration_workspace_free(w2);
}