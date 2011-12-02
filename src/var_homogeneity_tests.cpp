//#include <memory>
#include <list>
#include <cmath>
#include <iostream>
#include "var_homogeneity_tests.h"


//#include "bartlett_test.h"



//double var_homogneity_test(snp_var_data *, bool, char**);//invoke one ofvarianc homogeneity test (Bartlett's, Levene's and so on)
//double bartlett_test(snp_var_data * snp, bool all_genogroup_only);
//double bartlett_test(std::list<my_small_vector> * samples);

//double levene_test(snp_var_data * snp, bool all_genogroup_only);
//double levene_test(std::list<my_small_vector> * samples);




//
//___________________________________________________________
chisq_df var_homogneity_test(snp_var_data *snp, int  testname)
{
if(testname == bartlett) return bartlett_test(snp);
else if(testname == levene) return leven_test(snp);
else if(testname == likelihood) return likelihood_test(snp);
else if(testname == kolmogorov_smirnov) return kolmogorov_smirnov_test(snp);
else Rprintf("error: VarABEL: var_homogneity_test: Unexpected name of test. Only bartlett, levene, likelihood, kolmogorov_smirnov are alowed\n");
return chisq_df(NA_value, NA_value_int);
}
//___________________________________________________________


//___________________________________________________________
chisq_df var_homogneity_test(std::list<my_small_vector> * samples, short unsigned testname)
{
if(testname == bartlett) return bartlett_test(samples);
else if(testname == levene) return leven_test(samples);
else if(testname == likelihood) return likelihood_test(samples);
else if(testname == kolmogorov_smirnov) return kolmogorov_smirnov_test(samples);
else Rprintf("error: VarABEL: var_homogneity_test: Unexpected name of test. Only bartlett, levene, likelihood, kolmogorov_smirnov are alowed\n");
return chisq_df(NA_value, NA_value_int);
}
//___________________________________________________________






//_____________________________________________________________
// Bartlett test for homogeneity in variances
//http://en.wikipedia.org/wiki/Bartlett%27s_test
//retrun test chi2 or -1 in case of a problem
//___________________________________________________________
chisq_df bartlett_test(snp_var_data * snp)
{
//http://en.wikipedia.org/wiki/Bartlett%27s_test

static unsigned N; N=0;
static double Sp_2; Sp_2=0;
static double sum_ni_1lnSi2; sum_ni_1lnSi2=0;
static double sum_1__n_1; sum_1__n_1=0;
static double var;
static unsigned geno_group_num; geno_group_num=0;


for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
	if(snp->COUNTS[i] == 0 || is_na(snp->COUNTS[i])) continue;
	var = snp->SD[i]*snp->SD[i];
	
	N += snp->COUNTS[i];
	
	sum_ni_1lnSi2 += (snp->COUNTS[i]-1.)*log(var);
	sum_1__n_1 += 1./(snp->COUNTS[i]-1.);
	Sp_2 += (snp->COUNTS[i]-1.)*var;
	
	geno_group_num++;
	}

if(geno_group_num == 0 || geno_group_num == 1) return chisq_df(NA_value, NA_value_int);
Sp_2 /= N - geno_group_num; 

static double chisq;
static int df;
chisq = ((N - geno_group_num)*log(Sp_2) - sum_ni_1lnSi2)/(1 + (sum_1__n_1-1/(N-geno_group_num))/(3*(geno_group_num-1)) );
df = geno_group_num - 1;

return chisq_df(chisq, df);
}

//___________________________________________________________


//___________________________________________________________
chisq_df bartlett_test(std::list<my_small_vector> * samples)
	{
	static snp_var_data snp;
//	static double var;

	std::list<my_small_vector>::iterator it;
	
//	if(samples->size() != GENO_TYPES_NUM) 
//		{
//		Rprintf("warning: VarABEL: bartlett_test: Unexpected number of group\n");
//		return chisq_df(NA_value, NA_value);
//		}
	static unsigned i;
	i=0;
	for(it=samples->begin(); it!=samples->end(); ++it)
		{
		snp.SD[i] = my_var(&(*it));
		snp.COUNTS[i] = it->number;
		i++;
		}

	static chisq_df chi2;
	chi2 = bartlett_test(&snp);
	snp.reset();
	
	return chi2;

//	unsigned sample_amount = samples->size();
//
//	std::list<my_small_vector>::const_iterator it;
//
//	double Sp2=0;
//	double sum_ni_1lnSi2 = 0; //sum( (n_i-1)*ln(Si2) )
//	double sum_1__n_1=0;
//	double N=0;
//	double var_i=0;
//
//	for(it=samples->begin(); it!=samples->end(); ++it)
//		{
//		N += it->number;
//		var_i = var(*it);
//		sum_ni_1lnSi2 += (it->number - 1)*log(var_i);
//		sum_1__n_1 += 1./(it->number-1.);
//		Sp2 += (it->number - 1)*var_i;
//		} 
//
//	Sp2 = Sp2/(N-sample_amount);
//
//	double X2 = ((N - sample_amount)*log(Sp2) - sum_ni_1lnSi2)/(1 + (sum_1__n_1-1/(N-sample_amount))/(3*(sample_amount-1)) );
	}
//_____________________________________________________________













//___________________________________________________________
chisq_df leven_test(snp_var_data * snp)
{
//http://en.wikipedia.org/wiki/Levene%27s_test


static unsigned N; N=0; //total id's number
static VARIABLE_TYPE common_mean; common_mean=0;
static VARIABLE_TYPE var_sum; var_sum=0; //this is double sum in denumerator of the formula for Levene's test
static VARIABLE_TYPE var_mean; var_mean=0; //this is sum in numerator of the formula for Levene's test
static int geno_group_num; geno_group_num=0;




for(int i=0; i<GENO_TYPES_NUM ; i++)
	{
  if(snp->COUNTS[i] == 0 || is_na(snp->COUNTS[i])) continue;
	common_mean += snp->MEAN[i]*snp->COUNTS[i];
	N += snp->COUNTS[i];
	geno_group_num++;
	}


if(geno_group_num == 0 || geno_group_num==1) return chisq_df(NA_value, NA_value_int);




common_mean /= N;



//snp->SD[i]*snp->SD[i] - variance of trait's distribution in genotypic group group i



for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
  if(snp->COUNTS[i] == 0 || is_na(snp->COUNTS[i])) continue;
	var_sum  += snp->SD[i]*snp->SD[i]*snp->COUNTS[i];
	var_mean += snp->COUNTS[i]*(snp->MEAN[i] - common_mean)*(snp->MEAN[i] - common_mean);
	}



static double chisq;
static int df;

df = geno_group_num - 1;

chisq = (N-geno_group_num)*var_mean/(var_sum*df);


return chisq_df(chisq, df);

}




chisq_df leven_test(std::list<my_small_vector> * samples)
	{
	static snp_var_data snp;

	std::list<my_small_vector>::iterator it;

//	if(samples->size() != GENO_TYPES_NUM) 
//		{
//		Rprintf("warning: VarABEL: levene_test: Unexpected number of group (%i)\n", samples->size());
//		return chisq_df(NA_value, NA_value);;
//		}


	static unsigned i;
	i=0;
	for(it=samples->begin(); it!=samples->end(); ++it)
		{
		static double median;
//		median = my_mean(&(*it));
		median = my_median(&(*it));
		
		for(unsigned id=0 ; id<it->number ; id++)
			{
			it->vector[id] = fabs(it->vector[id] - median);
			}

		snp.MEAN[i] = my_mean(&(*it));
		snp.SD[i] = sqrt(my_var(&(*it)));
		snp.COUNTS[i] = it->number;
		i++;	
		}
	
	
	static chisq_df chi2;
	
	chi2 = leven_test(&snp);
	snp.reset();
	
	return chi2;
	}
//_____________________________________________________________





//_____________________________________________________________
// Likelihood test for homogeneity in variances
//retrun test chi2 or -1 in case of a problem

chisq_df likelihood_test(snp_var_data * snp)
	{
	return chisq_df(NA_value, NA_value_int);
	}

chisq_df likelihood_test(std::list<my_small_vector> * samples)
	{
	return chisq_df(NA_value, NA_value_int);
	}
//_____________________________________________________________


//_____________________________________________________________
// Kolmogorov test for homogeneity of variances
//retrun test chi2 or -1 in case of a problem

chisq_df kolmogorov_smirnov_test(snp_var_data * snp)
	{
  return chisq_df(NA_value, NA_value_int);	
	}

chisq_df kolmogorov_smirnov_test(std::list<my_small_vector> * samples)
	{
	return chisq_df(NA_value, NA_value_int);
	}
//_____________________________________________________________









