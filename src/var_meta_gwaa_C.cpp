//#=====================================================================================
//#
//#       Filename:  var_meta_gwaa_C.cpp
//#
//#    Description:  Function for meta analysis ov variance. Read plink format flat files and perform metanalysis.
//#
//#        Version:  1.0
//#        Created:  02-July-2009
//#       Revision:  none
//#
//#
//#         Author:  Maksim V. Struchalin
//#        Company:  ErasmusMC, Epidemiology, The Netherlands.
//#          Email:  m.struchalin@erasmusmc.nl
//#
//#=====================================================================================


//Description:
//
//Function
//"void var_meta_plink_C(const char** filenames, unsigned *file_amount_, const char **output_filename_, unsigned *skip_first_lines_amount_, char **delim_)"
//is wrtitten to meta analyze variance of snps. Input parameteres are 
//filenames - massive of file names where snp variances store. These files came from plink in cases of running it with key --qt-means like in bellow example
//plink --bfile $genotypes --allow-no-sex --pheno $phe_filename --all-pheno --missing-phenotype -999 --qt-means --out out --assoc
//input file format of such file is:
//CHR          SNP  VALUE      G11      G12      G22
//1   rs11497407   GENO      A/A      A/G      G/G
//1   rs11497407 COUNTS        0        2     5699
//1   rs11497407   FREQ        0 0.0003508   0.9996
//1   rs11497407   MEAN       NA  -0.2434 8.542e-05
//1   rs11497407     SD       NA   0.6545        1
//1   rs12565286   GENO      C/C      C/G      G/G
//1   rs12565286 COUNTS        3      459     5239
//1   rs12565286   FREQ 0.0005262  0.08051    0.919
//1   rs12565286   MEAN   0.2318 -0.02798 0.002319
//1   rs12565286     SD   0.1208    0.966    1.003
//
//Function var_meta_plink_C is robust enougph to little deviation from this format: 
//1) order of columns can be arbitrary
//2) there can be any other collumns. unnecessary columns will be ignored
//3) same is with unnecessary rowa. it is ignored if it is unnecessary


//input parameter "file_amount_" is number of files
//output_filename_ - output filename. All metaanalysed snps will be there
//skip_first_lines_amount_ - how many snps should be skipped in input files. 0 by default. In cases somebody wants use differing from plink format
//delim_ - which delim is used, blank space (' ') by default.





#include <map>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
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


#define VARIABLE_TYPE double
#define GENO_TYPES_NUM 3


#include <cmath>
#include "var_homogeneity_tests.h"

#include "supplementary_functions.h"
		
extern "C" {

#include "inverse_variance_metaanalysis.h"















//Start of function body 
//___________________________________________________________
//Function get names of files where are means and variance for each geno group. Return metanalysed variances and means.
//Function get as many functions as you have. Metanalysis goes file by file. 
//Files must have columns CHR, SNP, VALUE, G11, G12, G22. And each snp has to have parameters
//GENO, COUNTS, MEAN, SD (order isn't important). And don't argue with me!!!

void var_meta_gwaa_C(const char** filenames, unsigned *file_amount_, const char **output_filename_, unsigned *skip_first_lines_amount_, char **delim_,
											double* lambdas, unsigned* lambdas_NA,
											bool* exclude_whole_snp_if_number_less_threshold_,
											unsigned * threshold_,
											bool *do_all_warnings_output_,
											char* testname)
{
/*
bool exclude_whole_snp_if_number_less_threshold = *exclude_whole_snp_if_number_less_threshold_;
unsigned threshold = *threshold_;

bool do_all_warnings_output = *do_all_warnings_output_;

std::string output_filename = output_filename_[0];

std::ofstream warnings_file;
std::string warnings_filename = output_filename + ".warnings";
warnings_file.open(warnings_filename.c_str());
if(!warnings_file.is_open()){error("Can not open file %s\n", warnings_filename.c_str());}

		
unsigned file_amount = *file_amount_;
unsigned skip_first_lines_amount = *skip_first_lines_amount_;
char delim = *delim_[0];

double median; //inflation factor
std::string lambda_str;


Rprintf("\nChecking format of files...\n");
if(!check_files_format(filenames, file_amount, skip_first_lines_amount, delim))
	{
	Rprintf("Stop. Wrong file format.\n");
	for(unsigned i=0 ; i<=file_amount ; i++)
		{
		lambdas_NA[i]=1;
		}
	return;
	}
	
Rprintf("\nFiles are checked and OK.\nStart reading...\n\n");

		
Snp_store_type snps_data; //all analysed data will be here


//Rprintf("Each file has to contain columns with names %s, %s, %s, %s, %s, %s\n", chromosome_column_name.c_str(),
//																																								snpname_column_name.c_str(),
//																																								value_column_name.c_str(),
//																																								g11_column_name.c_str(),
//																																								g12_column_name.c_str(),
//																																								g22_column_name.c_str());



//read the first file
std::ifstream file;

std::stringstream chromosome_stream, snpname_stream, g11_value_stream, g12_value_stream, g22_value_stream;
std::stringstream num_to_string;
std::string value;
std::string str_from_stream;

std::vector<double> variance_tests_vec;//, lambda_vec;





for(unsigned file_num=0 ; file_num < file_amount ; file_num++)
	{
	Rprintf("\nProcessing file \"%s\"...\n", filenames[file_num]);
	
	warnings_file<<"Processing file \""<<filenames[file_num]<<"\"\n";
	warnings_file<<"___________________________________________\n";
	
	file.open(filenames[file_num]);
	if(!file.is_open()){error("Can not open file %s\n", filenames[file_num]);}

	//skip first line
	for(unsigned i=0 ; i<skip_first_lines_amount ; i++) 
		{
		getline(file, str_from_stream);
		if(file.eof()) {error("Tried to skip %i lines in file %s but there is %i at all ", skip_first_lines_amount, filenames[file_num], i);}
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

	if(CHR_position==-1) {error("Can not find column \"%s\"\n", chromosome_column_name.c_str());}
	if(SNP_position==-1) {error("Can not find column \"%s\"\n", snpname_column_name.c_str());}
	if(VALUE_position==-1) {error("Can not find column \"%s\"\n", value_column_name.c_str());}
	if(G11_position==-1) {error("Can not find column \"%s\"\n", g11_column_name.c_str());}
	if(G12_position==-1) {error("Can not find column \"%s\"\n", g12_column_name.c_str());}
	if(G22_position==-1) {error("Can not find column \"%s\"\n", g22_column_name.c_str());}



	//start reading data from file
	line_stream.str("");
	line_stream.clear(); //set eof to false
	
	static bool GENO_bool, COUNTS_bool, MEAN_bool, SD_bool;
	GENO_bool = COUNTS_bool = MEAN_bool = SD_bool = false;

	static unsigned excluded_snp_num;
	excluded_snp_num=0;

	static unsigned snp_number;
	snp_number=0;
	
	static unsigned current_line_number;
	current_line_number=0;

	static unsigned step;
	step = 100000;

	while(1) // one iteration = one snp
		{
		snp_var_data *snp = new snp_var_data;	

		snp_number++;
	
		while(1) //run this cycle untill a snp is ready
			{
			getline(file, str_from_stream); //get next line from file
			
			current_line_number++;
					
			if(file.eof()) break; //exit if EOF
				
			line_stream.str(str_from_stream);	
			
			for(int col=0; !line_stream.eof() ; col++ )  //start spliting current line
				{
				getline(line_stream, str_from_stream, delim);
			  if(str_from_stream.size() == 0) {col--; continue;}
				
				if(str_from_stream == "NA" || str_from_stream == "na" || str_from_stream == "nan" || str_from_stream == "NaN") 
					{
					num_to_string << NA_value;
					str_from_stream = num_to_string.str();
					num_to_string.str("");
					num_to_string.clear();
					}
			
				
				if(col == CHR_position)				 {chromosome_stream << str_from_stream;}
				else if(col == SNP_position) 	 {snpname_stream << str_from_stream;}
				else if(col == VALUE_position) {value = str_from_stream;}
				else if(col == G11_position)   {g11_value_stream << str_from_stream;}
				else if(col == G12_position)   {g12_value_stream << str_from_stream;}
				else if(col == G22_position)   {g22_value_stream << str_from_stream;}
				}
			//Now the line has been splited. Let's try to recognize what we splited and put it into snp info storage	

			
			//check wether we read same snp as in prevoius iteration. If this is new SNP then skip snp with incomplete info
			if(snp->snpname == "")
				{
				snp->snpname = snpname_stream.str();
				chromosome_stream >> snp->chromosome;
				}
			else if(snp->snpname != snpname_stream.str() && !(GENO_bool && COUNTS_bool && MEAN_bool && SD_bool))
				{
				//Rprintf("warning: Can not find complete information for SNP \"%s\". Line %i in file. SNP skiped\n", snp->snpname.c_str(), current_line_number);
				warnings_file<<"warning: Can not find complete information for SNP \""<<snp->snpname<<"\". Line "<<current_line_number<<" in file. SNP skiped\n";
				snp_number--;
				//Start reading new snp
				snp->reset();
				snp->snpname = snpname_stream.str();
				chromosome_stream >> snp->chromosome;
				GENO_bool=false, COUNTS_bool=false, MEAN_bool=false, SD_bool=false;
				}
			
			
			if(value == geno_value_name) {g11_value_stream >> snp->GENO[0]; g12_value_stream >> snp->GENO[1]; g22_value_stream >> snp->GENO[2]; GENO_bool = true;}
			else if(value == counts_value_name) {g11_value_stream >> snp->COUNTS[0]; g12_value_stream >> snp->COUNTS[1]; g22_value_stream >> snp->COUNTS[2]; COUNTS_bool = true;}
			else if(value == mean_value_name) {g11_value_stream >> snp->MEAN[0]; g12_value_stream >> snp->MEAN[1]; g22_value_stream >> snp->MEAN[2]; MEAN_bool = true;}
			else if(value == sd_value_name) {g11_value_stream >> snp->SD[0]; g12_value_stream >> snp->SD[1]; g22_value_stream >> snp->SD[2]; SD_bool = true;}

			g11_value_stream.str(""); g11_value_stream.clear();
			g12_value_stream.str(""); g12_value_stream.clear();
			g22_value_stream.str(""); g22_value_stream.clear();
			chromosome_stream.str(""); chromosome_stream.clear();
			snpname_stream.str(""); snpname_stream.clear();
			
			line_stream.str(""); line_stream.clear(); //set eof to false
			
			value="";


			

			if(GENO_bool && COUNTS_bool && MEAN_bool && SD_bool)
 				{
					
				if(!snp_filter(snp, warnings_file, exclude_whole_snp_if_number_less_threshold, threshold, do_all_warnings_output)) {excluded_snp_num++; break;}
				
				if(!unify_snp(snp, warnings_file)) {excluded_snp_num++; break;} //Start reading new SNP
				
				if(include_snp(&snps_data, snp, warnings_file, testname)) // put new snp into the storage. If this snp is there already tham metaanalyse it
					{	
					static chisq_df chisq;	
					chisq = var_homogneity_test(snp, testname);
					if(chisq.df != 1) {variance_tests_vec.push_back(chisq.chisq);}
					}
				else
					{
					excluded_snp_num++;
					}
				GENO_bool=false, COUNTS_bool=false, MEAN_bool=false, SD_bool=false;
				break; //Start reading new SNP
				}
		
			}


		if(snp_number % step == 0) 
			{
			Rprintf("%i SNPs done\n", snp_number);
			if(snp_number >= step*5) step *= 5;
			}


		if(file.eof()) {snp_number--; break;} //exit if EOF
		}
		
	file.close();
	file.clear();
	

	median = my_median(&variance_tests_vec);
	lambdas[file_num] = (is_na(median)? NA_value : median/median_table_chi2_df2);

	variance_tests_vec.clear();

	Rprintf("All SNPs done. Total amount of SNPs is %i. Excluded %i.\n", snp_number, excluded_snp_num);

	warnings_file<<"___________________________________________\n";
	warnings_file<<"End of file \""<<filenames[file_num]<<"\"\n\n";
	
	lambda_str = is_na(lambdas[file_num])? "NA" : double_2_str(lambdas[file_num]);
	Rprintf("Inflation factor for 2df variance homogeneity tests for snps from \"%s\" is lambda=%s \n", filenames[file_num], lambda_str.c_str());

	} // all files are read and snp pooled




Snp_store_type::iterator iter_map;
variance_tests_vec.clear();


//___________________________________________________________________________________
//perform variance homogeneity test for pooled data
	for(Snp_store_type::const_iterator i=snps_data.begin() ; i!=snps_data.end() ; ++i)
		{
		static chisq_df chisq;
		chisq = var_homogneity_test(i->second, testname);
		
		if(chisq.df == 1) i->second->Z = chisq.chisq;
		if(chisq.df == 2) i->second->Z_2df = chisq.chisq;
		
		if(!is_na(i->second->Z_2df)) variance_tests_vec.push_back(i->second->Z_2df);
		}
//___________________________________________________________________________________




//Obtain inflation factor
//___________________________________________________________________________________
median = my_median(&variance_tests_vec);
lambdas[file_amount] = (is_na(median)? NA_value : median/median_table_chi2_df2);
lambda_str = (is_na(lambdas[file_amount])? "NA":double_2_str(lambdas[file_amount]));
Rprintf("inflation factor for 2df variance homogeneity tests for pooled snps is lambda=%s \n", lambda_str.c_str());
//___________________________________________________________________________________







//get NA positions
//___________________________________________________________________________________
for(unsigned i=0 ; i<=file_amount ; i++)
	{
	if(is_na(lambdas[i])) {lambdas_NA[i]=1;}
	else {lambdas_NA[i]=0;}
	}
//___________________________________________________________________________________






//Save pooled means/variances and test statistics into files with extensions ".means" and ".tests" correspondingly
//___________________________________________________________________________________
std::string output_means_filename = output_filename + ".means";
std::string output_tests_filename = output_filename + ".tests";

Rprintf("\nwriting pooled means and variances into file %s\n", output_means_filename.c_str());
save_snps_data_into_file(&snps_data, output_means_filename.c_str(), delim);


Rprintf("\nwriting variance homogeneity tests results into file %s\n", output_tests_filename.c_str());
save_snps_tests_into_file(&snps_data, output_tests_filename.c_str(), delim);
//___________________________________________________________________________________




warnings_file.close();


Rprintf("Done\n");
*/
}
//End of var_meta_plink_C
//___________________________________________________________________________________













}
