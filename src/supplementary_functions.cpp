//#=====================================================================================
//#
//#       Filename:  supplementary_functions.cpp
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


#include "supplementary_functions.h"
#include "inverse_variance_metaanalysis.h"
#include <stdlib.h>
#include <math.h>







//___________________________________________________________
//Set genotypes in common view GA -> AG
bool unify_snp(snp_var_data* snp, std::ofstream & warnings_file)
{


static std::string alleles_snp;
alleles_snp="";

for(unsigned i=0 ; i<GENO_TYPES_NUM ; i++) {alleles_snp += snp->GENO[i];}//collect all alleles 

alleles_snp = get_uniq_symbols(alleles_snp);
for(unsigned i=0 ; i<alleles_snp.size() ; i++) if(alleles_snp[i] == '/') {alleles_snp.erase(i,1);}

if(alleles_snp.size() > 2 && alleles_snp.size() <=1)
	{
	//Rprintf("SNP %s has more than 2 alleles (it has %s). SNP skiped.\n", snp->snpname.c_str(), alleles_snp.c_str());
	warnings_file<<"SNP "<<snp->snpname<<" has more than 2 alleles (it has "<<alleles_snp<<"). SNP skiped.\n";
	return false;	
	}

if(alleles_snp.size() <=1)
	{
	//Rprintf("SNP %s has less than 2 alleles (it has %s). SNP skiped.\n", snp->snpname.c_str(), alleles_snp.c_str());
	warnings_file<<"SNP "<<snp->snpname<<" has less than 2 alleles (it has "<<alleles_snp<<"). SNP skiped.\n";
	return false;	
	}


if(alleles_snp[0]=='0' && alleles_snp[1]=='0')
	{
	//Rprintf("SNP %s has undefind allels. SNP skiped..\n", snp->snpname.c_str());
	warnings_file<<"SNP "<<snp->snpname<<" has undefind allels. SNP skiped\n";
	return false;
	}




//exclude those geno group which has only one id
for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{

	if(is_na(snp->SD[i]) || is_na(snp->MEAN[i]) || is_na(snp->COUNTS[i]) || snp->COUNTS[i]<=1 )
		{
		snp->SD[i] = NA_value;
		snp->MEAN[i] = NA_value;
		snp->COUNTS[i] = 0;
		continue;
		}

	static VARIABLE_TYPE var;
	var = snp->SD[i]*snp->SD[i];	
	if(var <= 1.E-32)
		{
		snp->SD[i] = NA_value;
		snp->MEAN[i] = NA_value;
		snp->COUNTS[i] = 0;
		
		//Rprintf("warning: genotypic group %s in snp %s has too small variance (variance=%f). This genotypic group is excluded from analysis.\n", 
		warnings_file<<"warning: genotypic group "<<snp->GENO[i]<<" in snp "<<snp->snpname<<" has too small variance (variance="<<var<<"). This genotypic group is excluded from analysis.\n"; 
		
		}
	}






std::sort(alleles_snp.begin(), alleles_snp.end());


if(alleles_snp[0] == '0')
	{
	alleles_snp[0] = alleles_snp[1];
	alleles_snp[1] = '0';
	}

static std::string geno[3];


geno[0].push_back(alleles_snp[0]);
geno[0].push_back('/');
geno[0].push_back(alleles_snp[0]);
		
geno[1].push_back(alleles_snp[0]);
geno[1].push_back('/');
geno[1].push_back(alleles_snp[1]);

geno[2].push_back(alleles_snp[1]);
geno[2].push_back('/');
geno[2].push_back(alleles_snp[1]);


static std::string geno_hemozyg_switched;
geno_hemozyg_switched.push_back(geno[1][2]);
geno_hemozyg_switched.push_back('/');
geno_hemozyg_switched.push_back(geno[1][0]);


if(snp->GENO[0] == geno[0] && snp->GENO[1] == geno[1] && snp->GENO[2] == geno[2])
	{
	geno_hemozyg_switched = "";
	geno[0] = "";
	geno[1] = "";
	geno[2] = "";
	return true; //nothing to change
	}


static snp_var_data snp_new;

//snp_new->snpname = snp->snpname;
//snp_new->chromosome = snp->chromosome;


for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
	if(snp->GENO[i] == geno[0])
		{snp_new.GENO[0] = geno[0]; snp_new.COUNTS[0]=snp->COUNTS[i]; snp_new.MEAN[0]=snp->MEAN[i]; snp_new.SD[0]=snp->SD[i];}
	if(snp->GENO[i] == geno[1] || snp->GENO[i]==geno_hemozyg_switched) 
		{snp_new.GENO[1] = geno[1]; snp_new.COUNTS[1]=snp->COUNTS[i]; snp_new.MEAN[1]=snp->MEAN[i]; snp_new.SD[1]=snp->SD[i];}
	if(snp->GENO[i] == geno[2]) 
		{snp_new.GENO[2] = geno[2]; snp_new.COUNTS[2]=snp->COUNTS[i]; snp_new.MEAN[2]=snp->MEAN[i]; snp_new.SD[2]=snp->SD[i];}
	}

geno_hemozyg_switched = "";
geno[0] = "";
geno[1] = "";
geno[2] = "";









for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
	snp->GENO[i] = snp_new.GENO[i];
	snp->COUNTS[i] = snp_new.COUNTS[i];
	snp->MEAN[i] = snp_new.MEAN[i];
	snp->SD[i] = snp_new.SD[i];
	}

snp_new.reset();




return true;
}
//___________________________________________________________




//___________________________________________________________
//put new snp into the storage. If this snp is there already tham metaanalyse it
bool include_snp(Snp_store_type * snps_storage, snp_var_data* snp, std::ofstream & warnings_file, char *testname)
{


Snp_store_type::iterator iter_map = snps_storage->find(snp->snpname);

//char delim=' ';



if(iter_map == snps_storage->end()) 
	{
	snps_storage->insert(std::pair<std::string, snp_var_data*>(snp->snpname, snp));
	}
else
	{
	snp_var_data* snp_current = iter_map->second;
	static bool is_snp_ok;
	is_snp_ok = check(snp_current, snp, warnings_file);
	
	if(!is_snp_ok) {return false;}
	
	(*snps_storage)[snp->snpname.c_str()] = snp_var_meta(iter_map->second, snp, testname);
	
	//bless god souls of these objects...
	delete snp_current;
	delete snp;
	}

return true;
}
//___________________________________________________________



//metaanalysis of data from two snps
//___________________________________________________________
snp_var_data* snp_var_meta(snp_var_data* snp1, snp_var_data* snp2, char* testname)
{
snp_var_data* snp_meta = new snp_var_data;



snp_meta->chromosome = snp1->chromosome;
snp_meta->snpname = snp1->snpname;

VARIABLE_TYPE meta_mean_beta[GENO_TYPES_NUM], meta_mean_se[GENO_TYPES_NUM];
	

static unsigned one = 1;

for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
	//metanalysis for variance

	static std::string meta_codding;
			
	snp_meta->GENO[i] = snp1->GENO[i];
	
	if(snp1->COUNTS[i] != 0	&& snp2->COUNTS[i] != 0 )
		{

		static VARIABLE_TYPE SD_of_the_mean_snp1, SD_of_the_mean_snp2;
		
		SD_of_the_mean_snp1 = snp1->SD[i]/sqrt(double(snp1->COUNTS[i]));
		SD_of_the_mean_snp2 = snp2->SD[i]/sqrt(double(snp2->COUNTS[i]));

		//metanalysis for mean
		inverse_variance_metaanalysis(&snp1->MEAN[i], &snp2->MEAN[i],
						 &SD_of_the_mean_snp1, &SD_of_the_mean_snp2,
							&one,
							meta_mean_beta,
							meta_mean_se, testname);
	
		snp_meta->MEAN[i] = meta_mean_beta[0];
		snp_meta->COUNTS[i] = snp1->COUNTS[i] + snp2->COUNTS[i];
		snp_meta->SD[i] = meta_mean_se[0]*sqrt(double(snp_meta->COUNTS[i]));
		}
	else
		{
		if(snp1->COUNTS[i] == 0)
			{
			snp_meta->SD[i] = snp2->SD[i];
			snp_meta->COUNTS[i] = snp2->COUNTS[i];
			snp_meta->MEAN[i] = snp2->MEAN[i];
			}
		else if(snp2->COUNTS[i] == 0)
			{
			snp_meta->SD[i] = snp1->SD[i];
			snp_meta->COUNTS[i] = snp1->COUNTS[i];
			snp_meta->MEAN[i] = snp1->MEAN[i];
			}
		else
			{
			error("Upss.. :-) something strange occured... wrong file format probably. Create small example of your data where you have an error and send it to developer.\n");
			}
		}
	}
return snp_meta;		
}
//___________________________________________________________







//___________________________________________________________
bool is_na(const VARIABLE_TYPE val, const VARIABLE_TYPE na_reference)
{
static VARIABLE_TYPE delta = fabs(na_reference/1E6);
if(val > na_reference-delta && val < na_reference+delta) return true;
else return false;
}
//___________________________________________________________





bool is_na(const int val, const int na_reference) //is numerical value recognized as NA 
{
if(val == na_reference) return true;
else return false;
}



//save all data into plink like format
//___________________________________________________________
void save_snps_data_into_file(Snp_store_type *snps_data, const char *output_filename, char delim)
{
//print header
	
const unsigned pp_maxsnp=10;
		
std::ofstream file;
	

file.open(output_filename);


if(!file.is_open()){error("Can not open file %s\n", output_filename);}



//file<<setiosflags(std::ios::right);

file.precision(precision_output);

file << std::setw(4) << chromosome_column_name.c_str() << delim
		 << std::setw(pp_maxsnp) << snpname_column_name.c_str() << delim
		 << std::setw(6) << value_column_name.c_str() << delim
		 << std::setw(8) << g11_column_name.c_str() << delim
		 << std::setw(8) << g12_column_name.c_str() << delim
		 << std::setw(8) << g22_column_name.c_str() << "\n";


		
static Snp_store_type::iterator iter_map;

for(Snp_store_type::const_iterator i=snps_data->begin() ; i!=snps_data->end() ; ++i)
	{
	//geno:
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << geno_value_name << delim
			 << std::setw(8) << i->second->GENO[0] << delim
			 << std::setw(8) << i->second->GENO[1] << delim
			 << std::setw(8) << i->second->GENO[2] << "\n";



	//counts:
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << counts_value_name << delim;

		 if(is_na(i->second->COUNTS[0]))
			 {
			 file <<  std::setw(8) << "NA" << delim;
			 }
			else
			 {
			 file <<  std::setw(8) << std::setiosflags(std::ios_base::dec) << i->second->COUNTS[0] << delim;
			 }

		 if(is_na(i->second->COUNTS[1]))
			 {
			 file <<  std::setw(8) << "NA" << delim;
			 }
			else
			 {
			 file <<  std::setw(8) << std::setiosflags(std::ios_base::dec) << i->second->COUNTS[1] << delim;
			 }

		 if(is_na(i->second->COUNTS[2]))
			 {
			 file <<  std::setw(8) << "NA" << delim;
			 }
			else
			 {
			 file <<  std::setw(8) << std::setiosflags(std::ios_base::dec) << i->second->COUNTS[2] << '\n';
			 }



	//freq:
	static VARIABLE_TYPE total_id_num;
 	total_id_num = i->second->COUNTS[0] + i->second->COUNTS[1] + i->second->COUNTS[2];	
	
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << freq_value_name << delim
			 << std::setw(8) << (is_na(i->second->COUNTS[0])? "NA":double_2_str(i->second->COUNTS[0]/total_id_num)) << delim
			 << std::setw(8) << (is_na(i->second->COUNTS[1])? "NA":double_2_str(i->second->COUNTS[1]/total_id_num)) << delim
			 << std::setw(8) << (is_na(i->second->COUNTS[2])? "NA":double_2_str(i->second->COUNTS[2]/total_id_num)) << "\n";


	//mean:
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << mean_value_name << delim
			 << std::setw(8) << (is_na(i->second->MEAN[0])? "NA":double_2_str(i->second->MEAN[0])) << delim
			 << std::setw(8) << (is_na(i->second->MEAN[1])? "NA":double_2_str(i->second->MEAN[1])) << delim
			 << std::setw(8) << (is_na(i->second->MEAN[2])? "NA":double_2_str(i->second->MEAN[2])) << "\n";

	//sd:
	file << std::setw(4) << i->second->chromosome << delim
		   << std::setw(pp_maxsnp) << i->second->snpname << delim 
			 << std::setw(6) << sd_value_name << delim
			 << std::setw(8) << (is_na(i->second->SD[0])? "NA":double_2_str(i->second->SD[0])) << delim
			 << std::setw(8) << (is_na(i->second->SD[1])? "NA":double_2_str(i->second->SD[1])) << delim
			 << std::setw(8) << (is_na(i->second->SD[2])? "NA":double_2_str(i->second->SD[2])) << "\n";


	}


file.close();
}
//___________________________________________________________




//___________________________________________________________
void save_snps_tests_into_file(Snp_store_type *snps_data, const char *output_filename, char delim)
{
//print header
	
//const unsigned pp_maxsnp=10;
		
std::ofstream file;
	

file.open(output_filename);


if(!file.is_open()){error("Can not open file %s\n", output_filename);}



//file<<setiosflags(std::ios::right);

file.precision(precision_output);

file << chromosome_column_name << snpname_column_name << delim << "Z" << delim << "Z_2df\n";


		
static Snp_store_type::iterator iter_map;

for(Snp_store_type::const_iterator i=snps_data->begin() ; i!=snps_data->end() ; ++i)
	{
	//geno:
	file << i->second->chromosome << delim
			 << i->second->snpname << delim 
			 << (is_na(i->second->Z)? "NA":double_2_str(i->second->Z)) << delim
			 << (is_na(i->second->Z_2df)? "NA":double_2_str(i->second->Z_2df)) << "\n"; 
	}


file.close();
}
//___________________________________________________________



//convert from double to string
//___________________________________________________________
std::string double_2_str(double val, const unsigned precision)
{
static std::stringstream stream;
stream.str(""); stream.clear();
stream.precision(precision);
stream << val;
return stream.str();
}
//___________________________________________________________



//___________________________________________________________
//Check whether snps have same genotypes and in same columns. If one of snp has allele 0 but another have real than replace 0 by real one.
bool check(snp_var_data* snp2, snp_var_data* snp1, std::ofstream & warnings_file)
{
//check snp name
if(snp1->snpname != snp2->snpname) 
	{
	error("snp_var_meta: unexpected error; atempt to pool two different snps");
	}

//check chromosome name
if(snp1->chromosome != snp2->chromosome)
	{
	//Rprintf("warning: SNP %s has different chromosome number in different files. Previos value is %i, current one is %i. Value %i is used.\n",
	//			 	snp2->snpname.c_str(), snp2->chromosome, snp1->chromosome, snp2->chromosome);
	
	warnings_file<<"warning: SNP "<<snp2->snpname<<" has different chromosome number in different files. Previos value is "<<snp2->chromosome<<
			", current one is "<<snp1->chromosome<<". Value "<<snp2->chromosome<<" is used.\n";
	}



//check coddings
if(snp1->GENO[1] == snp2->GENO[1])
	{
	return true;
	}



//snp1 - snp from next cohort. If genotypes from snp1 and snp2 don't match than skip snp1

if(snp1->GENO[1][0] == snp2->GENO[1][0]) //A1==A2?
	{
	if(snp1->GENO[1][2] == '0') //B1==0?
		{
		if(snp2->GENO[1][2] != '0') //B2!=0?
			{
			//situation when snp1 has A/A, A/0, 0/0 and snp2 has A/A, A/B, B/B. Replace 0 by B in snp1
			snp1->GENO[1][2] = snp2->GENO[1][2];
			snp1->GENO[2][0] = snp2->GENO[1][2];
			snp1->GENO[2][2] = snp2->GENO[1][2];
			return true;
			}
		}
	else //B1!=0!!!
		{
		if(snp2->GENO[1][2] == '0') // B2==0?
			{
			//situation when snp1 has A/A, A/B, B/B and snp2 has A/A, A/0, 0/0. Replace 0 by B in snp1
			snp2->GENO[1][2] = snp1->GENO[1][2];
			snp2->GENO[2][0] = snp1->GENO[1][2];
			snp2->GENO[2][2] = snp1->GENO[1][2];
			return true;
			}
		else
			{
			//B1!=0 and B2!=0 therefore codings defer. Skip snp1
			//Rprintf("warning: snp %s has different genotypes in current and previous cohort. Current one is %s, previos - %s. The current one is skiped.\n", 
			//				snp1->snpname.c_str(), snp1->GENO[1].c_str(), snp2->GENO[1].c_str());
			
			warnings_file<<"warning: snp "<<snp1->snpname<<
					" has different genotypes in current and previous cohort. Current one is "<<snp1->GENO[1]<<", previous - "<<snp2->GENO[1]<<
					". The current one is skiped.\n";
			return false;
			}
		}
	
	}
else if(snp1->GENO[1][0] == snp2->GENO[1][2])
		{
		//situation when snp1 has B/B, B/C, C/C and snp2 has A/A, A/B, B/B. C can be 0. Let's check it
		if(snp1->GENO[1][2] == '0')
			{
			//situation when snp1 has B/B, B/0, 0/0 and snp2 has A/A, A/B, B/B.
			snp1->GENO[1][2] = snp2->GENO[1][0];
			snp1->GENO[2][0] = snp2->GENO[1][0];
			snp1->GENO[2][2] = snp2->GENO[1][0];
			unify_snp(snp1, warnings_file);
			return true;
			}
		else //B1!=0!!!
			{
			//it could be snp1 has B/B, B/A, A/A and snp2 has A/A, A/0, 0/0. But snp1 has been sorted already => second allele of snp2 is C!=0
			//Rprintf("warning: snp %s has different genotypes in current and previous cohort. current is %s, previos is %s. The current one is skiped.\n", 
			//				snp1->snpname.c_str(), snp1->GENO[1].c_str(), snp2->GENO[1].c_str());
			
			warnings_file<<"warning: snp "<<snp1->snpname<<
					" has different genotypes in current and previous cohort. Current one is "<<snp1->GENO[1]<<", previous is "<<snp2->GENO[1]<<
					". The current one is skiped.\n"; 
					
			return false;
			}


		}
else if(snp1->GENO[1][2] == snp2->GENO[1][0])
	{
	//situation when snp1 has A/A, A/B, B/B and snp2 has B/B, B/0, 0/0. Replace 0 by B in snp1
	snp2->GENO[1][2] = snp1->GENO[1][0];
	snp2->GENO[2][0] = snp1->GENO[1][0];
	snp2->GENO[2][2] = snp1->GENO[1][0];
	unify_snp(snp2, warnings_file);
	
	}
else
	{
	//A1!=A2 and A1!=B2
	//Rprintf("warning: snp %s has different genotypes in current and previos cohort. current is %s, previos is %s. The current one is skiped.\n", 
	//				snp1->snpname.c_str(), snp1->GENO[1].c_str(), snp2->GENO[1].c_str());
	warnings_file<<"warning: snp "<<snp1->snpname
							 <<" has different genotypes in current and previous cohort. Current one is "<<snp1->GENO[1]
							 <<", previos is "<<snp2->GENO[1]
							 <<". The current one is skiped.\n"; 
	
	return false;
	}

return true;
}
//___________________________________________________________



//___________________________________________________________
// input parameter is "A/AA/GG/G", output AG/
std::string get_uniq_symbols(std::string alleles_snp)
{
std::string uniqe_symbols="";
int size = alleles_snp.size();

for(int i=0 ; i<size ; i++)
	{
	static char symbol;
 	symbol = alleles_snp[i];
	
	static int size_uniqe;
	size_uniqe = uniqe_symbols.size();

	static bool flag;
	flag=false;
	for(int j=0 ; j<size_uniqe ; j++)
		{
		if(uniqe_symbols[j] == symbol) {flag = true; break;}
		}
	
	if(!flag) {uniqe_symbols += symbol;}

	}


return uniqe_symbols;
}
//___________________________________________________________




//___________________________________________________________
double my_median(std::vector<double> * vec)
{
unsigned size = vec->size();

if(size == 0) return NA_value;
if(size == 1) return (*vec)[0];

//for(int i=0 ; i<vec->size() ; i++)
//	{
//	}

std::sort(vec->begin(), vec->end());

//for(int i=0 ; i<vec->size() ; i++)
//	{
//	}


static double median;

if(size % 2 < 1E-12) {median =  ((*vec)[size/2-1] + (*vec)[size/2])/2;} //odd number
else {median = (*vec)[(size-1)/2];} //even number


return median;

}
//_____________________________________________________________



//___________________________________________________________
double my_median(my_small_vector * vec)
{
unsigned size = vec->number;

if(size == 0) return NA_value;
if(size == 1) return vec->vector[0];


qsort(vec->vector, vec->number, sizeof(double), compare_doubles);


static double median;

if(size % 2 < 1E-12) { median = (vec->vector[size/2-1] + vec->vector[size/2])/2;} //odd number
else {median = vec->vector[(size-1)/2];} //even number


return median;
}
//_____________________________________________________________

//_____________________________________________________________
double my_var(my_small_vector * vec)
	{
	double sum = 0;
	double mean = my_mean(vec);

	if(vec->number <2) error("error: var: Phenotypic data contains less than two subjects.");


	for(unsigned i=0; i<vec->number ; i++)
		{
		sum += (vec->vector[i]-mean)*(vec->vector[i]-mean);
		}

	return sum/(vec->number-1);
	}

//_____________________________________________________________



//_____________________________________________________________
// returun mean of the array or exit(1) in case of problem
double my_mean(my_small_vector * vec)
	{
	static double mean;
  mean = 0;
	

	if(vec->number == 0) {Rprintf("error: get_mean: sample does not have any element"); return NA_value;}

	for(unsigned i=0; i<vec->number ; i++)
		{
		mean += vec->vector[i];
		}

	return mean/double(vec->number);
	}
//_____________________________________________________________


//___________________________________________________________


//___________________________________________________________
bool snp_filter(snp_var_data* snp,  std::ofstream & warnings_file, bool exclude_whole_snp, unsigned threshold, bool do_warnings_output)
{


for(int i=0 ; i<GENO_TYPES_NUM ; i++)
	{
  if(snp->COUNTS[i]<=0) continue;
	if(unsigned(snp->COUNTS[i]) < threshold)
		if(exclude_whole_snp)
			{
			if(do_warnings_output)
				{
				warnings_file<<"warning: SNP \""<<snp->snpname<<"\" has been excluded because of number of ids in genotyoic group "<<snp->GENO[i]<<" less than "<<threshold<<'\n';
				}
			return false;
			}
		else
			{
			if(do_warnings_output)
				{
				warnings_file<<"warning: genotypic group "<<snp->GENO[i]<<" in SNP \""<<snp->snpname<<"\" has been excluded because of number of ids less than "<<threshold<<'\n';
				}
			snp->SD[i] = NA_value;
			snp->MEAN[i] = NA_value;
			snp->COUNTS[i] = 0;
			}	
	
	}

return true;

}
//___________________________________________________________








//___________________________________________________________
bool check_files_format(const char** filenames, unsigned file_amount, unsigned skip_first_lines_amount, char delim)
{

std::stringstream chromosome_stream, snpname_stream, g11_value_stream, g12_value_stream, g22_value_stream;
std::ifstream file;
std::string str_from_stream;

std::stringstream num_to_string;

for(unsigned file_num=0 ; file_num < file_amount ; file_num++)
	{
	Rprintf("\nChecking file \"%s\"...\n", filenames[file_num]);
	
	
	file.open(filenames[file_num]);
	if(!file.is_open()){Rprintf("Can not open file %s\n", filenames[file_num]); return false;}

	//skip first line
	for(unsigned i=0 ; i<skip_first_lines_amount ; i++) 
		{
		getline(file, str_from_stream);
		if(file.eof()) {Rprintf("Tried to skip %i lines in file %s but there is %i at all ", skip_first_lines_amount, filenames[file_num], i);return false;}
		}

	//read header and determine position of our columns

	int CHR_position=-1, VALUE_position=-1, SNP_position=-1, G11_position=-1, G12_position=-1, G22_position=-1; //0 means te first column

	getline(file, str_from_stream);
	std::stringstream line_stream(str_from_stream);

	if(file.eof()) break;

	for(unsigned col=0 ; !line_stream.eof() ; col++ )
		{
		getline(line_stream, str_from_stream, delim);
		if(str_from_stream.size() == 0) {col--; continue;}
	

		if(str_from_stream == chromosome_column_name) {CHR_position=col;}
		else if(str_from_stream == snpname_column_name) {SNP_position=col;}
		else if(str_from_stream == value_column_name) {VALUE_position=col;}
		else if(str_from_stream == g11_column_name) {G11_position=col;}
		else if(str_from_stream == g12_column_name) {G12_position=col;}
		else if(str_from_stream == g22_column_name) {G22_position=col;}
		}

	if(CHR_position==-1) {Rprintf("Can not find column \"%s\".\n", chromosome_column_name.c_str()); return false;}
	if(SNP_position==-1) {Rprintf("Can not find column \"%s\".\n", snpname_column_name.c_str()); return false;}
	if(VALUE_position==-1) {Rprintf("Can not find column \"%s\".\n", value_column_name.c_str()); return false;}
	if(G11_position==-1) {Rprintf("Can not find column \"%s\".\n", g11_column_name.c_str()); return false;}
	if(G12_position==-1) {Rprintf("Can not find column \"%s\".\n", g12_column_name.c_str()); return false;}
	if(G22_position==-1) {Rprintf("Can not find column \"%s\".\n", g22_column_name.c_str()); return false;}


	Rprintf("File \"%s\" is OK\n", filenames[file_num]);
		
	file.close();
	file.clear();
	} // all files are checked


return true;
}





//The function break_trait_up_into_groups takes phenotype and snp, break the phenotype up on three genotypic groups.
//_________________________________________________________________________________________________
void break_trait_up_into_groups(std::list<my_small_vector> *trait_groups, double* snp, double *trait, long unsigned* nids, int analys_type, int * is_trait_na)
{

std::vector<double> NA, AA, AB, BB; // here we will store trait for different genotype group

	for(unsigned id=0 ; id < *nids ; id++)
		{
		if(is_trait_na[id] == 1) continue;

	//get bestguess:
	//________________		
	static int snp_value;
	if(snp[id]>=0. && snp[id]<0.5 ) snp_value=0;
	if(snp[id]>=0.5 && snp[id]<=1.5 ) snp_value=1;
	if(snp[id]>1.5 && snp[id]<=2. ) snp_value=2;
	//________________		


		//spread ids trait among genotype group
		switch(snp_value)
			{
//			case 0:
//				{
//				//NA.push_back(trait[id]);
//				}
			break;
			case 0:
				{
				AA.push_back(trait[id]);
				}
			break;
			case 1:
				{
				AB.push_back(trait[id]);
				}
			break;
			case 2:
				{
				BB.push_back(trait[id]);
				}
			break;
			default:
				{
				Rprintf("error: VarABEL: Unexpected genotype code has been detected (%i, %i). Only 0, 1, 2, NA are alowed\n", &snp[id], snp[id]);
				return;	
				}
			break;
			}
		}

	
	
//_________________________________________________
//1df conversion:
	//b.insert(b.end(), a.begin(), a.end());
  //enum {AAvsABvsBB, AAvsABandBB, ABvsAAandBB, BBvsAAandAB};

	if(analys_type == AAvsABandBB)
		{
		AB.insert(AB.end(), BB.begin(), BB.end());
		BB.clear();
		}
	else if(analys_type == ABvsAAandBB)
		{
		AA.insert(AA.end(), BB.begin(), BB.end());
		BB.clear();
		}
	else if(analys_type == BBvsAAandAB)
		{
		AB.insert(AB.end(), AA.begin(), AA.end());
		AA.clear();
		}






//_________________________________________________

//	unsigned NA_size = NA.size();
	unsigned AA_size = AA.size();
	unsigned AB_size = AB.size();
	unsigned BB_size = BB.size();

	double /**na,*/ *aa, *ab, *bb;


		if(AA_size > 1)
			{
			aa = new double[AA_size];
			for(unsigned i=0 ; i<AA_size ; i++)
				{
				aa[i] = AA[i];
				}
			trait_groups->push_back(my_small_vector(aa, AA_size));
			}



		if(AB_size > 1)
			{
			ab = new double[AB_size];
			for(unsigned i=0 ; i<AB_size ; i++)
				{
				ab[i] = AB[i];
				}
			trait_groups->push_back(my_small_vector(ab, AB_size));
			}
	

	
		if(BB_size > 1)
			{
			bb = new double[BB_size];
			for(unsigned i=0 ; i<BB_size ; i++)
				{
				bb[i] = BB[i];
				}
			trait_groups->push_back(my_small_vector(bb, BB_size));
			}

		
//if(AA_size > 1) delete[] aa;
//if(AB_size > 1) delete[] ab;
//if(BB_size > 1) delete[] bb;

//AA.clear();
//AB.clear();
//BB.clear();
		
		
		
		
return;
}//end of function break_trait_up_into_groups	
//_________________________________________________________________________________________________



//________________________________________________________
// For qsort function
int compare_doubles(const void *a, const void *b)
{
double* arg1 = (double*) a;
double* arg2 = (double*) b;
if( *arg1 < *arg2 ) return -1;
else if( *arg1 == *arg2 ) return 0;
else return 1;
}
//________________________________________________________




