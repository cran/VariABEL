//#=====================================================================================
//#
//#       Filename:  gtps_container.cpp
//#
//#    Description:  Container for storaging genotype data. Use in function for merging of two snp.data class objects.
//#
//#        Version:  1.0
//#        Created:  18-March-2008
//#       Revision:  none
//#       
//#
//#         Author:  Maksim V. Struchalin, Yurii S. Aulchenko
//#        Company:  ErasmusMC, Epidemiology & Biostatistics Department, The Netherlands.
//#          Email:  m.struchalin@@erasmusmc.nl, i.aoultchenko@erasmusmc.nl
//#
//#=====================================================================================

#include "gtps_container.h"
#include <cstdlib>  
#include <math.h>




gtps_container::gtps_container(char * gtps_array_, char * strand_array_, char * coding_array_, unsigned id_numbers_, unsigned snp_numbers_)
{
do_we_have_strand_and_codding_arrays = true; //of course
rearrangement_array = new unsigned[4];
rearrangement_array[0] = 6;
rearrangement_array[1] = 4;
rearrangement_array[2] = 2;  
rearrangement_array[3] = 0;





gtps_array = gtps_array_;	
strand_array = strand_array_;
coding_array = coding_array_;
id_numbers = id_numbers_;
snp_numbers = snp_numbers_;
nbytes_for_one_snp = unsigned(ceil(double(id_numbers)/4.) + 0.5);


}
//------------------------------------------------------------------



//------------------------------------------------------------------
gtps_container::gtps_container(char * gtps_array_, unsigned id_numbers_, unsigned snp_numbers_)
{
do_we_have_strand_and_codding_arrays = true; //of course
rearrangement_array = new unsigned[4];
rearrangement_array[0] = 6;
rearrangement_array[1] = 4;
rearrangement_array[2] = 2;  
rearrangement_array[3] = 0;





gtps_array = gtps_array_;	
id_numbers = id_numbers_;
snp_numbers = snp_numbers_;
nbytes_for_one_snp = unsigned(ceil(double(id_numbers)/4.) + 0.5);
}
//------------------------------------------------------------------






gtps_container::~gtps_container(void)
{
delete[] rearrangement_array;
}




//------------------------------------------------------------------
void gtps_container::get_our_byte_number_and_local_person_number(unsigned id_position, unsigned snp_position)
{
our_byte_number = int(ceil(id_position/4.)+0.5) + (snp_position-1)*nbytes_for_one_snp;   //What byte is our id in? 

local_person_number = id_position - ((our_byte_number-(snp_position-1)*nbytes_for_one_snp-1)*4); //What position of our id is in byte?
}
//------------------------------------------------------------------









//------------------------------------------------------------------
char gtps_container::get(unsigned id_position, unsigned snp_position)
{
get_our_byte_number_and_local_person_number(id_position, snp_position); //calculate our_byte_number value

//char our_byte_vallue = gtps_array[our_byte_number-1];



return char((gtps_array[our_byte_number-1] >> rearrangement_array[local_person_number-1]) & 3); 
}
//------------------------------------------------------------------




//------------------------------------------------------------------
//Attention! This class does not carry about allocated in this function memory.
//In order to avoid memory leak "delete" must be performed for every array created in this function.
char* gtps_container::get_gtps_array_for_snp(unsigned snp_position)
{
char* gtps_for_one_snp = new char(nbytes_for_one_snp);
get_our_byte_number_and_local_person_number(1, snp_position); //calculate our_byte_number value



for(unsigned int i=0 ; i<nbytes_for_one_snp; i++)
	{
	gtps_for_one_snp[i]=gtps_array[our_byte_number-1+i];
	}


return gtps_for_one_snp;
}
//------------------------------------------------------------------









char gtps_container::get_strand(unsigned snp_position)
{
if(do_we_have_strand_and_codding_arrays) return strand_array[snp_position-1]; 
else Rprintf("gtps_container::get_strand: You can not get strand since you create object with constructor gtps_container(char * gtps_array_raw, unsigned id_numbers, unsigned snp_numbers)");
return char(0);
}



char gtps_container::get_coding(unsigned snp_position)
{
if(do_we_have_strand_and_codding_arrays) return coding_array[snp_position-1]; 
else Rprintf("gtps_container::get_strand: You can not get strand since you create object with constructor gtps_container(char * gtps_array_raw, unsigned id_numbers, unsigned snp_numbers)");
return char(0);
}


//------------------------------------------------------------------
void gtps_container::set(unsigned id_position, unsigned snp_position, char data)
{

		
//static const char clear_info_for_person[]={63, 207, 243, 252};
//static const char clear_info_for_person[]={'&#255;', '&#207;', '&#243;', '&#252;'};
static const char clear_info_for_person[]={char(63), char(207), char(243), char(252)};

get_our_byte_number_and_local_person_number(id_position, snp_position);



//char tmp = (gtps_array[our_byte_number]&clear_info_for_person[local_person_number-1]) | data;




gtps_array[our_byte_number-1] = (gtps_array[our_byte_number-1]&clear_info_for_person[local_person_number-1]) | (data&3)<<rearrangement_array[local_person_number-1];
}
//------------------------------------------------------------------




/*
gtps_container::checking(unsigned id_position, unsigned snp_position)
{
if (id_position > id_numbers || id_position<0)
	{
	}

if(snp_position > snp_position || snp_position <0)
	{
	}
}
*/



