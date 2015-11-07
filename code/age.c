#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"

/**@file age.c Returns time to present for a given redshift.
 *
 * @brief For a given redshift, returns the time from that redshift until
 * redshift zero:
 * \f$H_0t_0=\int_0^z\frac{dz}{(1+z)\sqrt{(1+z)^2(1+z\Omega_m)-z(2+x)\Omega_{\Lambda}}}\f$
 *
 * Returns Age in code units/Hubble_h*/

double time_to_present(double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double time, result, abserr;
 
  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;

  /* I think that the 1.0/Hubble here should be any small number, such
     as 1e-6 - it is the aboslute error accuracy required.  PAT */
  gsl_integration_qag(&F, 1.0 / (z + 1), 1.0, 1.0 / Hubble,
                      1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  time = 1 / Hubble * result;

  gsl_integration_workspace_free(workspace);

  
  
  return time;
}


double integrand_time_to_present(double a, void *param)
{
  return 1 / sqrt(Omega / a + (1 - Omega - OmegaLambda) + OmegaLambda * a * a);
}
