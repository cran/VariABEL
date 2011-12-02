//#=====================================================================================
//#
//#       Filename:  var_homogeneity_test_C.cpp
//#
//#    Description:  Function to perform variance homogeneity test in R. 
//#
//#        Version:  0.1
//#        Created:  27-Apr-2010
//#       Revision:  none
//#				 Modifed:  
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl
//#
//#=====================================================================================


#include "var_homogeneity_tests.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <Rinternals.h>
#include <fstream>
#include <algorithm>

#include <vector>

#include <R.h>  // to include Rconfig.h

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("pkg", String)
// replace pkg as appropriate
#else
#define _(String) (String)
#endif



#include <Rinternals.h>
#include <Rdefines.h>

#include <cmath>
#include <list>

#include "gtps_container.h"
#include "supplementary_functions.h"
#include "constants.h"

#include <libintl.h>

#include "linear_regression.h"






extern "C" {

//Tne function which are invoked by DatABEL iterator:

//Variance homogeneity test. The function is invoked from R. The function only prepare trait data and invoke one of the available variance homogeneity tests.
//_________________________________________________________________________________________________


void variance_homogeneity_test_C(double* snp, double *trait, double *design_matrix_copy, double *design_matrix, int *p_, long unsigned* nids_, double *betas, double *se, double * chi2, int * df, double * residuals, int *analys_type_, int * is_trait_na, int *testname, /*auxiliary variables:*/ double *qty, int *jpvt, double *qraux, double *work, double *v, double *x_for_ch2inv)
{
		
std::list<my_small_vector> trait_groups;
static chisq_df chisq;



if(*testname != svlm)
	{
	break_trait_up_into_groups(&trait_groups, snp, trait, nids_, *analys_type_, is_trait_na); 

	//Start analysis of homogeneity
	chisq = var_homogneity_test(&trait_groups, *testname);
	*chi2 = chisq.chisq;
	*df = chisq.df;
	}
else
	{
	int p;
	p = *p_;
	static int ii=0;
	
	ii++;

	long int nids = *nids_;



	//fillfull the last covariates by the values from snp
	for(int i=nids*(p-1) ; i<nids*p ; i++) 
		{
	
		design_matrix_copy[i]=snp[i-nids*(p-1)];
		}
	
//	for(int i=0 ; i<nids*p ; i++) 
//		{
//		}
		

//	for(int i=0 ; i<nids ; i++)
//		{
//		}





	//drop out the NA values
	//___________________________________________
//	long unsigned int nids_nona=0;
//	for(int i=0 ; i<nids ; i++)
//		{
//		nids_nona += is_trait_na[i];
//		}
//	nids_nona = nids - nids_nona;	

	long unsigned int nids_nona=0;
	for(int column=0 ; column<p ; column++)
		{
		for(int id=0; id<nids ; id++)
			{
			if(is_trait_na[id] == 0) 
				{
				design_matrix[nids_nona] = design_matrix_copy[id + column*nids];
				
				nids_nona++;	
				}
			}	
		}

	nids_nona=0;
	for(int id=0; id<nids ; id++)
		{
		if(is_trait_na[id] == 0)
			{
			trait[nids_nona] = trait[id];
			nids_nona++;
			}
		}

	

	//___________________________________________

	
	//Get rid of major SNP effect
	//___________________________________________
	
	int p_disp = 2;
	for(long unsigned int i=0; i<nids_nona ; i++) design_matrix_copy[i] = design_matrix[i];
	for(long unsigned int i=0; i<nids_nona ; i++) design_matrix_copy[i+nids_nona] = design_matrix[i +(p-1)*nids_nona];

//	for(int i=0; i<2*nids_nona ; i++)
//		{
//		}

//	for(int i=0; i<nids_nona ; i++)
//		{
//		}


 
 linear_regression(
				        /*input variables:*/     trait, design_matrix, p_, &nids_nona,
								/*return variables:*/    betas, se, residuals,
								/*auxiliary variables:*/ qty, jpvt, qraux, work, v, x_for_ch2inv);

	//___________________________________________

//  double beta_snp   = betas[1];
//	double sebeta_snp = se[1];
	
		




	//___________________________________________

 	double betas_disp[2], se_disp[2];

	for(int i=0; i<nids ; i++) {residuals[i] = residuals[i]*residuals[i];}
	linear_regression(
				        /*input variables:*/     residuals, design_matrix_copy, &p_disp, &nids_nona,
								/*return variables:*/    betas_disp, se_disp, residuals,
								/*auxiliary variables:*/ qty, jpvt, qraux, work, v, x_for_ch2inv);

	//___________________________________________




	*chi2 = (betas_disp[1]/se_disp[1])*(betas_disp[1]/se_disp[1]);
	*df = 1;
	betas[p] = betas_disp[1];
  se[p]	   = se_disp[1];
	
	}

return;
} // end of function bartlett_test
//_________________________________________________________________________________________________







} // end of extern "C"


//extern "C" {
//
//
//SEXP variance_homogeneity_test_C_old_data_type(SEXP set_, SEXP trait_, SEXP nids_, SEXP nsnps_, SEXP analys_type_, SEXP testname_)
//{
//char *set = (char*)(RAW(set_));
//double *trait = REAL(trait_);
//unsigned nids = INTEGER_VALUE(nids_);
//int nsnps = INTEGER_VALUE(nsnps_);
//const char * analys_type_char = CHARACTER_VALUE(analys_type_);
//const char * testname_char = CHARACTER_VALUE(testname_);
//
//
//
//gtps_container Set(set, NULL, NULL, nids, nsnps); //creat object to facilitate working with set1
//
//int *snp = new int[nids];
//double chi2;
//double *chi2_vec;
//int df;
//int *df_vec;
//
//chi2_vec = new double[nsnps];
//df_vec = new int[nsnps];
//int* is_trait_na = new int[nids];
//
//
////Fill is_trait_na. 1 means missing vallue, 0 - non missing
////_________________________________________________________
//for(unsigned id_counter=0 ; id_counter<nids ; id_counter++)
//	{
//	if(ISNAN(trait[id_counter])) 
//		{
//		is_trait_na[id_counter]=1;
//		}
//	else
//		{
//		is_trait_na[id_counter]=0;
//		}
//	}
////_________________________________________________________
//
//
////_________________________________________________________
//analysis_type_enum analys_type;
//
//if(strcmp(analys_type_char, "AAvsABvsBB") == 0)
//	{
//	analys_type = AAvsABvsBB;
//	}
//else if(strcmp(analys_type_char, "AAvsABandBB") == 0)
//	{
//	analys_type = AAvsABandBB;
//	}
//else if(strcmp(analys_type_char, "ABvsAAandBB") == 0)
//	{
//	analys_type = ABvsAAandBB;
//	}
//else if(strcmp(analys_type_char, "BBvsAAandAB") == 0)
//	{
//	analys_type = BBvsAAandAB;
//	}
//else
//	{
//	error("Unknown name of analysis type. Available are AAvsABvsBB, AAvsABandBB, ABvsAAandBB, BBvsAAandAB");
//	}
////_________________________________________________________
//
//
//
////_________________________________________________________
//testname_enum testname;
//if(strcmp(testname_char, "bartlett") == 0)
//	{
//	testname = bartlett;
//	}
//else if(strcmp(testname_char, "levene") == 0)
//	{
//	testname = levene;
//	}
//else if(strcmp(testname_char, "likelihood") == 0)
//	{
//	testname = likelihood;
//	}
//else if(strcmp(testname_char, "kolmogorov_smirnov") == 0)
//	{
//	testname = kolmogorov_smirnov;
//	}
//else
//	{
//	error("Unknown name of test name. Available are bartlett, levene, likelihood, kolmogorov_smirnov");
//	}
////_________________________________________________________
//
//unsigned step=10000;
//
//
//for(int snp_position=1 ; snp_position<=nsnps ; snp_position++)
//	{
//	for(int id_counter=1 ; id_counter<=nids ; id_counter++) {snp[id_counter-1] = int(Set.get(id_counter, snp_position));}
//	variance_homogeneity_test_C(snp, trait, &nids, &chi2, &df, analys_type, is_trait_na, testname);
//	chi2_vec[snp_position-1] = chi2;
//	df_vec[snp_position-1] = df;
//	
//	if(snp_position % step == 0)
//		{
//		Rprintf("%d SNPs done\n", snp_position);
//		if(snp_position >= step*5) step *= 5;
//		}
//	}
//
////void variance_homogeneity_test_C(int* snp, double *trait, unsigned* nids, double * chi2, double * df, char* analys_type_, int * is_trait_na, char* testname)
//
//
//
//SEXP chi2_df_results_R_LIST;
//PROTECT(chi2_df_results_R_LIST = NEW_LIST(2)); 
//
//SET_VECTOR_ELT(chi2_df_results_R_LIST, 0, NEW_NUMERIC(nsnps));
//SET_VECTOR_ELT(chi2_df_results_R_LIST, 1, NEW_INTEGER(nsnps));
//
//for(unsigned i = 0; i < nsnps; i++)
//	{
//	NUMERIC_POINTER(VECTOR_ELT(chi2_df_results_R_LIST, 0))[i] = chi2_vec[i];
//	}
//
//for(unsigned i = 0; i < nsnps; i++)
//	{
//	INTEGER_POINTER(VECTOR_ELT(chi2_df_results_R_LIST, 1))[i] = df_vec[i];
//	}
//
//UNPROTECT(1);
//
////SEXP chi2_results_R;
////PROTECT(chi2_results_R = allocVector(REALSXP, 2*nsnps));
////double *chi2_results_R_double = REAL(chi2_results_R);
////
////for(unsigned i = 0; i < nsnps; i++)
////	{
////	chi2_results_R_double[i] = chi2_vec[i];
////	}
////
////for(unsigned i = nsnps; i < 2*nsnps; i++)
////	{
////	chi2_results_R_double[i] = df_vec[i];
////	}
////UNPROTECT(1);
//
//
//delete snp;
//delete chi2_vec;
//delete df_vec;
//delete is_trait_na;
//
//return(chi2_df_results_R_LIST);
//}
//
//
//
//} //end of extern "C"
