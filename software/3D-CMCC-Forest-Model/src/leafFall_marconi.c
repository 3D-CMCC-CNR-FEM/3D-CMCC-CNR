#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpfit.h"
#include "types.h"
#include "constants.h"

struct vars_struct {
	double *x;
	double *y;
	double *ey;
};


// Gaussian fitting of VpdSat data on the function of the type : y = p1 * exp(-(x-p2)/p3)^2 (Gaussian Function
// we are interested in fitting the VPD sat data to find out the value of the p2 parameter value, used to asses the first day of senescence on the
// basis of the VPDsat peak value assumption
// The residuals are defined as (y[i] - f(x[i]))/y_error[i]
// Original public domain version by B. Garbow, K. Hillstrom, J. More'
// (Argonne National Laboratory, MINPACK project, March 1980)
// Tranlation to C Language by S. Moshier (moshier.net)
// Enhancements, documentation and packaging by C. Markwardt
// (comparable to IDL fitting routine MPFIT
// see http://cow.physics.wisc.edu/~craigm/idl/idl.html)
// Copyright (C) 2003, 2004, 2006, 2007, 2009, 2010 Craig B. Markwardt

void get_vpsat(CELL * c,  int day, int month, int years, YOS *yos, int i)
{
	MET_DATA *met;
	// check parameters
	met = (MET_DATA*) yos[years].m;

	c->vpSat[i] = 0.35*(0.6107*exp(17.38*(met[month].d[day].tavg)/(239+(met[month].d[day].tavg))))+0.154;
	//soil_Log("\nvpSat: %g", c->vpSat[i]);
}

int gaussian_x_y(int m, int n, double *p, double *dy,
		double **dvec, void *vars)
{
	int i;
	struct vars_struct *v = (struct vars_struct *) vars;
	double *x, *y, *ey;
	double xc;
	/* Retrieve values of x, y and y_error from private structure */

	x = v->x;
	y = v->y;
	ey = v->ey;


	for (i=0; i<m; i++) {
		xc = x[i]-p[1];
		dy[i] = (y[i] - p[0] * exp(-pow((xc/p[2]),2)))/ey[i];
		//dy[i] = (y[i] - p[1]*exp(-0.5*xc*xc/sig2) - p[0])/ey[i];
	}


	return 0;
}

int leaffalMarconi(CELL * c)
{
	/* Test harness routine, which contains test data, invokes mpfit() */
	/* X - independent variable */
	double x[365];
	/* Y - measured value of dependent quantity */
	double y[365];
	double ey[365];   /* Measurement uncertainty - initialized below */

	//	double x[] = {-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
	//			-2.8394297E-01,1.3416969E+00,1.3757038E+00,
	//			-1.3703436E+00,4.2581975E-02,-1.4970151E-01,
	//			8.2065094E-01};
	//	double y[] = {-4.4494256E-02,8.7324673E-01,7.4443483E-01,
	//			4.7631559E+00,1.7187297E-01,1.1639182E-01,
	//			1.5646480E+00,5.2322268E+00,4.2543168E+00,
	//			6.2792623E-01};
	//	double ey[10];

	double p[] = {0.0, 1.0, 1.0, 1.0};       /* Initial conditions */
	double pactual[] = {0.0, 4.70, 0.0, 0.5};/* Actual values used to make data*/
	double perror[4];                   /* Returned parameter errors */
	mp_par pars[4];			   /* Parameter constraints */
	int i;
	struct vars_struct v;  /* Private data structure */
	int status;
	mp_result result;


	//filling the x array with DOY values
	for(i=0; i<365; i++)
	{
		x[i] = i+1;
	}

	//filling the y array with VPDsat values
	for(i = 0; i< 365; i++)
	{
		//use the formula to evaluate VPDsat using surface temperature; compute it elsewhere, without including here the nth for cycle
		y[i] = c->vpSat[i];
	}

	memset(&result,0,sizeof(result));       /* Zero results structure */
	result.xerror = perror;

	memset(pars,0,sizeof(pars));        /* Initialize constraint structure */
	/* No constraints */
	for (i=0; i<365; i++) ey[i] = 0.5;   /* Data errors */



	/* Fill private data structure */
	v.x = x;
	v.y = y;
	v.ey = ey;

	/* Call fitting function for 10 data points and 2 parameters */
	status = mpfit(gaussian_x_y, 365, 4, p, pars, 0, (void *) &v, &result);

	//soil_Log("\n*** testlinfit status = %d\n", status);
	/* ... print or use the results of the fitted parametres p[] here! ... */
	//Marconi: turn on the following function to print the parameters
	//printresult(p, pactual, &result);

	return 0;

}

/* Simple routine to print the fit results */
void printresult(double *x, double *xact, mp_result *result)
{
	int i;

	if ((x == 0) || (result == 0)) return;
	soil_Log("  CHI-SQUARE = %f    (%d DOF)\n",
			result->bestnorm, result->nfunc-result->nfree);
	soil_Log("        NPAR = %d\n", result->npar);
	soil_Log("       NFREE = %d\n", result->nfree);
	soil_Log("     NPEGGED = %d\n", result->npegged);
	soil_Log("     NITER = %d\n", result->niter);
	soil_Log("      NFEV = %d\n", result->nfev);
	soil_Log("\n");
	if (xact) {
		for (i=0; i<result->npar; i++) {
			soil_Log("  P[%d] = %f +/- %f     (ACTUAL %f)\n",
					i, x[i], result->xerror[i], xact[i]);
		}
	} else {
		for (i=0; i<result->npar; i++) {
			soil_Log("  P[%d] = %f +/- %f\n",
					i, x[i], result->xerror[i]);
		}
	}

}


