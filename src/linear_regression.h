#ifndef SMV_LINEAR_REGRESSION_H
#define SMV_LINEAR_REGRESSION_H



//#=====================================================================================
//#
//#       Filename:  linear_regression.h
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



extern "C" {

extern int dqrls_(double *x, int *n, int *p, double *y, int *ny, double *tol, double *b, double *rsd, double *qty, int *k, int *jpvt, double *qraux, double *work);
extern void mych2inv_(double *x, int *ldx, int *n, double *v, int *info);

}

void linear_regression(
				        /*input variables:*/     double *y, double *design_matrix, int *p_/*cov num + 1*/, long unsigned* num_observations,
								/*return variables:*/    double *betas, double *se, double *residuals,
								/*auxiliary variables:*/ double *qty, int *jpvt, double *qraux, double *work, double *v, double *x_for_ch2inv);

#endif
