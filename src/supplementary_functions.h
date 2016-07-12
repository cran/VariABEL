//#=====================================================================================
//#
//#       Filename:  supplementary_functions.h
//#
//#    Description:  Supplementary functions for variance homogeneity test and for metaanalysis.
//#
//#        Version:  1.0
//#        Created:  28-Apr-2010
//#       Revision:  none
//#
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl
//#
//#=====================================================================================

#ifndef SMV_SUPPLEMENTARY_FUNCTIONS_H
#define SMV_SUPPLEMENTARY_FUNCTIONS_H

#include <map>
#include <iostream>
#include <list>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <algorithm>

#include <vector>

#include <R.h>
#include <Rinternals.h>
#undef length




#define VARIABLE_TYPE double
#define GENO_TYPES_NUM 3



#include "constants.h"

const std::string chromosome_column_name = "CHR";
const std::string snpname_column_name = "SNP";
const std::string value_column_name = "VALUE";
const std::string g11_column_name = "G11";
const std::string g12_column_name = "G12";
const std::string g22_column_name = "G22";

const std::string geno_value_name = "GENO";
const std::string counts_value_name = "COUNTS";
const std::string freq_value_name = "FREQ";
const std::string mean_value_name = "MEAN";
const std::string sd_value_name = "SD";





struct snp_var_data;
typedef std::map<std::string, snp_var_data*> Snp_store_type; //use std::string as key because for some reason it doesn't want to work with const char*



int compare_doubles(const void *a, const void *b);




class my_small_vector
	{
	
	public:
		my_small_vector(double * vector_, unsigned long number_)
			{
			vector=vector_;
			number=number_;
			}

		my_small_vector(const my_small_vector& p)
			{
			number = p.number;
			
			vector = new double[number];
			for(int i=0 ; i<number ; i++)
				{
				vector[i] = p.vector[i];
				}
			}


		~my_small_vector(void)
			{
			if(vector != NULL) delete[] vector;
			}

	double * vector; 
	long number; //amount of cells in vector


	};

//all info for a SNP is here
//_________________________________________________________
struct snp_var_data
	{
	snp_var_data() //set all variables to zero
		{
		reset();
		}	
			
	inline void reset() 
		{
		for(int i=0 ; i<GENO_TYPES_NUM ; i++)
			{
			GENO[i] = "NA";
			COUNTS[i] = NA_value_int;
			MEAN[i] = NA_value;
			SD[i]=NA_value;
			Z=NA_value;
			}
		chromosome=NA_value_int;
		snpname="";
		}
	
	std::string snpname;	
	std::string GENO[GENO_TYPES_NUM];
	int COUNTS[GENO_TYPES_NUM];
	VARIABLE_TYPE MEAN[GENO_TYPES_NUM];
	VARIABLE_TYPE SD[GENO_TYPES_NUM];
	int chromosome;
	double Z; //homogeneity test
	double Z_2df; //homogeneity test 2df only
	};
//_________________________________________________________



bool include_snp(Snp_store_type *, snp_var_data*, std::ofstream & warnings_file, char *testname); //include snp into common storage
snp_var_data* snp_var_meta(snp_var_data* , snp_var_data*, char* testname); //metaanalysis of two snps
bool is_na(const VARIABLE_TYPE val, const VARIABLE_TYPE na_reference=NA_value); //is numerical value recognized as NA 
bool is_na(const int val, const int na_reference=NA_value_int); //is numerical value recognized as NA 
void save_snps_data_into_file(Snp_store_type *snps_data, const char *output_filename, char delim); //save all snps into flat file in plink like format
void save_snps_tests_into_file(Snp_store_type *snps_data, const char *output_filename, char delim); //save all snps into flat file in plink like format
std::string double_2_str(VARIABLE_TYPE val, const unsigned precision=precision_output); // convert double to string
bool check(snp_var_data* snp1, snp_var_data* snp2, std::ofstream & warnings_file);

bool unify_snp(snp_var_data* snp, std::ofstream & warnings_file);
std::string get_uniq_symbols(std::string alleles_snp);

double my_median(my_small_vector * vec);
double my_median(std::vector<double> * vec);
double my_var(my_small_vector * vec);
double my_mean(my_small_vector * vec);

bool snp_filter(snp_var_data* snp,  std::ofstream & warnings_file, bool exclude_whole_snp, unsigned threshold, bool do_warnings_output);
bool check_files_format(const char** filenames, unsigned file_amount, unsigned skip_first_lines_amount, char delim);

void break_trait_up_into_groups(std::list<my_small_vector> *trait_groups, double* snp, double *trait, long unsigned* nids, int analys_type, int * is_trait_na);


const short unsigned int AAvsABvsBB = 0;
const short unsigned int AAvsABandBB = 1;
const short unsigned int ABvsAAandBB = 2;
const short unsigned int BBvsAAandAB = 3;

const short unsigned int bartlett = 0;
const short unsigned int levene = 1;
const short unsigned int likelihood = 2;
const short unsigned int kolmogorov_smirnov = 3;
const short unsigned int svlm=4;




#endif
