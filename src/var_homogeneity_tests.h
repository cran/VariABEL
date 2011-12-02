#ifndef SMV_VAR_HOMOGENEITY_TESTS_H
#define SMV_VAR_HOMOGENEITY_TESTS_H


#include <list>
#include "supplementary_functions.h"

#include "constants.h"

class chisq_df
	{
	public:
		chisq_df(const double chisq_, const int df_) 
			{
			chisq=chisq_;
			df=df_;
			}

		chisq_df() 
			{
			chisq=NA_value;
			df=NA_value_int;
			}
	
		chisq_df& operator=(const chisq_df& a) 
			{
			chisq = a.chisq;
			df    = a.df;
			return *this;
			}
	
		double chisq;
		int df;
	};
		



chisq_df var_homogneity_test(snp_var_data *, short unsigned int );//invoke one ofvarianc homogeneity test (Bartlett's, Levene's and so on)
chisq_df var_homogneity_test(std::list<my_small_vector> * samples, short unsigned int);

chisq_df bartlett_test(snp_var_data * snp);
chisq_df bartlett_test(std::list<my_small_vector> * samples);

chisq_df leven_test(snp_var_data * snp);
chisq_df leven_test(std::list<my_small_vector> * samples);

chisq_df likelihood_test(snp_var_data * snp);
chisq_df likelihood_test(std::list<my_small_vector> * samples);

chisq_df kolmogorov_smirnov_test(snp_var_data * snp);
chisq_df kolmogorov_smirnov_test(std::list<my_small_vector> * samples);

double get_mean(my_small_vector vec);
double var(my_small_vector vec);

extern "C" {
void variance_homogeneity_test_C(double* snp, double *trait, double *design_matrix_geno_means, double *design_matrix, int *p, long unsigned* nids, double *betas, double *se, double * chi2, int * df, double * residuals, int *analys_type_, int * is_trait_na, int *testname, /*auxiliary variables:*/ double *qty, int *jpvt, double *qraux, double *work, double *v, double *x_for_ch2inv);
}

#endif
