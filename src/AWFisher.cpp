/*******************************************************************************
pval is a matrix saved by row. pval = t(p.values)
weights is a matrix saved by row. The input should be all 0.

AW_weight() is for weight.matrix = TRUE.

g++ -O2 -fPIC -shared AW_weight.cpp pchisq.c -o AW_weight.so
********************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sort.h"
#include "pchisq.h"


#define        LOG_SQRT_PI     0.5723649429247000870717135 /* log (sqrt (pi)) */
#define        I_SQRT_PI       0.5641895835477562869480795 /* 1 / sqrt (pi) */
#define        BIGX            1000.0         /* max value to represent exp (x) */
#define CONST_SQRT2		(1.4142135623730951)

//#define        ex(x)             (((x) < -BIGX) ? 0.0 : exp (x))



__inline double _normcdf(double x)
{
	return .5 * (1. + erf(x / CONST_SQRT2));
}

double pchisq(double x, int df)
{
       double a, y, s;
       double e, c, z;
       int even;     /* true if df is an even number */

	   if (x <= 0.0 || df < 1)
	   {
		   return 1.;
	   }

       a = 0.5 * x;
	   even = !(df & 1);

	   //if (df > 1)
       y = exp (-a);
       s = (even ? y : (2.0 * _normcdf(-sqrt (x))));

	   if (df > 2)
	   {
		   x = 0.5 * (df - 1.0);
		   z = (even ? 1.0 : 0.5);

		   if (a > BIGX)
		   {
			   e = (even ? 0.0 : LOG_SQRT_PI);
			   c = log(a);
			   while (z <= x)
			   {
				   e = log(z) + e;
				   s += exp (c*z - a - e);
				   z += 1.0;
			   }
			   return s;
		   }
		   else
		   {
			   e = (even ? 1.0 : (I_SQRT_PI / sqrt(a)));
			   c = 0.0;
			   while (z <= x)
			   {
				   e = e * (a / z);
				   c = c + e;
				   z += 1.0;
			   }
			   return (c * y + s);
		   }
	   }
	   else
	   {
		   return s;
	   }
}


extern "C" void AWpvalue(double *best_stat, int *sum_weight, int *weights, 
	double *pval, int *nrow, int *ncol)
{
	int nr, nc, i, j, sw;
	double *pv, *ptwice_pcumlog, twice_pcumlog, stat_new, best;
	int *o, *pw;

	nr = *nrow;
	nc = *ncol;
	//nc1 = nc - 1;

	o = (int *)malloc(sizeof(int) * nc);
	ptwice_pcumlog = (double *)malloc(sizeof(double) * nc);

	pv = pval;
	pw = weights;
	for (i = 0; i < nr; ++i)
	{
		QuickSort1(nc, pv, o);

		sw = 0;

		twice_pcumlog = -2. * log(pv[o[0]]);
		best = pchisq(twice_pcumlog, 2);
		for (j = 1; j < nc; ++j)
		{
			twice_pcumlog -= 2. * log(pv[o[j]]);
			stat_new = pchisq(twice_pcumlog, j + j + 2);
			if (stat_new < best)
			{
				best = stat_new;
				sw = j;
			}

		}
		best_stat[i] = -log(best);
		sum_weight[i] = sw;
		
		for (j = 0; j <= sw; ++j)
		{
			pw[o[j]] = 1;
		}

		pv += nc;
		pw += nc;
	}

	free(ptwice_pcumlog);
	free(o);
}
