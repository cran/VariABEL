//#=====================================================================================
//#
//#       Filename:  inverse_variance_metaanalysis.h
//#
//#    Description:  Function for meta analysis with inverse variance method
//#
//#        Version:  1.0
//#        Created:  06-July-2009
//#       Revision:  none
//#				Modifed:  26-Apr-2010
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl
//#
//#=====================================================================================

#ifndef SMV_DOMETA_C_H
#define SMV_DOMETA_C_H


extern "C" {
void inverse_variance_metaanalysis(double *beta_set1, double *beta_set2,
              double *sebeta_set1, double *sebeta_set2,
              unsigned *num,
              double *mbeta,
              double *mse,
							char *testname);


}
#endif

