//#=====================================================================================
//#
//#       Filename:  inverse_variance_metaanalysis.cpp
//#
//#    Description:  Function for meta analysis with inverse variance method. 
//#
//#        Version:  1.0
//#        Created:  06-July-2009
//#       Revision:  none
//#				 Modifed:  26-Apr-2010
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl
//#
//#=====================================================================================




#include<iostream>
//#include<math>

#include <Rinternals.h>

#include <R.h>  // to include Rconfig.h 

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("pkg", String)
// replace pkg as appropriate 
#else
#define _(String) (String)
#endif

extern "C" {
#include "inverse_variance_metaanalysis.h"


void inverse_variance_metaanalysis(double *beta_set1, double *beta_set2,
			 				double *sebeta_set1, double *sebeta_set2,
			 				unsigned *num, //number of betas for current cohort
							double *mbeta,
							double *mse,
							char *testname)
{

unsigned num_el = *num;
double se_set1, se_set2, wt2_set1, wt2_set2, invsumwt2;



	for(unsigned i=0 ; i<num_el ; i++)
		{
		se_set1 = sqrt(sebeta_set1[i]*sebeta_set1[i]);
		se_set2 = sqrt(sebeta_set2[i]*sebeta_set2[i]);
	
		wt2_set1 = 1./(sebeta_set1[i]*sebeta_set1[i]);
		wt2_set2 = 1./(sebeta_set2[i]*sebeta_set2[i]);

		invsumwt2 = 1./(wt2_set1+wt2_set2);
		mbeta[i] = (beta_set1[i]*wt2_set1 + beta_set2[i]*wt2_set2)*invsumwt2;
		mse[i] = sqrt(invsumwt2);
		}

}// end of function 



}//end of extern "C"
