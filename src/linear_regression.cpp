#include "linear_regression.h"
#include <math.h>
#include "iostream"


//#=====================================================================================
//#
//#       Filename:  linear_regression.cpp
//#
//#    Description:  The function "linear_regression" computes least squares solutions to the system (1). Based on "dqrfit" and "ch2inv" fortran subroutine. 
//#									 (1) design_matrix * betas = y
//#
//#        Version:  0.1
//#        Created:  2-Dec-2010
//#       Revision:  none
//#				 Modifed:  
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl, m.v.struchalin@mail.ru
//#
//#=====================================================================================



//Manual of input parameters:
// "y" - dependable variable - array with length "num_observations"
// "design_matrix" - design (module) matrix with dimension ("num_observations" x "p"). The first column is filled with one's (for intercept). The following columns are covariates.
// "p" - number of covariates plus 1. "p" is number of columns in the design matrix.
// "num_observations" - number of observation. "num_observations" is length of "y" and number of rows of design matrix.
// "betas" - solution of equation (1). The array of length "p" which contains effects of intercect and covariates.
// "se" - solution of equation (1). The array of length "p" which contains standart errors of effects of intercect and covariates.
// "residuals" - residuals of y ("residuals"="y" - "betas"*"design_matrix"). 
//
// auxiliary variables: explanation of "qty", "jpvt", "qraux", "work" can be found in manual for dqrls (http://svn.r-project.org/R/trunk/src/appl/dqrls.f)
// "qty" - array of length "num_observations"
// "jpvt" - array of length "p"
// "qraux" - array of length "p"
// "work" - array of length "2*p"
// "v" - array of length "2*p". 
// "x_for_ch2inv" - array of length "2*p"


void linear_regression(
				/*input variables:*/     double *y, double *design_matrix, int *p_/*cov num + 1*/, long unsigned* num_observations, 
				/*return variables:*/    double *betas, double *se, double *residuals,
				/*auxiliary variables:*/ double *qty, int *jpvt, double *qraux, double *work, double *v, double *x_for_ch2inv)
{
	int n=*num_observations;
	int p = *p_;
//	int i_start = n*(p-1);
//	int i_stop  = n*p;

	static int ny = 1; //numeber of dependable variables (y)
	static double tol=1e-07;
	static int k;




	dqrls_(design_matrix, &n, &p, y, &ny, &tol, betas, residuals, qty, &k, jpvt, qraux, work);

	for(int col=0 ; col<p ; col++)
		{
		for(int row=0 ; row<p ; row++)
			{
			x_for_ch2inv[row + p*col] = design_matrix[row+n*col];
			}
		
		}

	int info;
	mych2inv_(x_for_ch2inv, &p, &p, v, &info);


	//calculate rss (residuals square sum)
	double rss=0;
	for(int i=0 ; i<n ; i++)
		{
		rss += residuals[i]*residuals[i];
		}
	
	static double rdf;
	rdf = n - p;
	static double resvar;
 	resvar = rss/rdf;


	for(int i=0 ; i<p ; i++)
		{
		se[i] = sqrt(v[i*p+i] * resvar); //diagonal ellements of v multiply by resvar
		}


}


