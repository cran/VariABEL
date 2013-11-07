#include "Rstuff.h"
#include "iterator_functions.h"
#include "var_homogeneity_tests.h"
#include <R.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif

	// For testing purposes, tried to move some functions+wrappers to GenABEL: failed
	// Product function + wrapper
	double prod(double *mydata, unsigned int size) {
		double prodtotal = mydata[0];
		for (unsigned int i = 1; i < size; i++) {
			prodtotal *= mydata[i];
		}
		return prodtotal;
	}
	void prodWrapper(double *indata, unsigned long int indataHeight,
			unsigned long int indataWidth,
			double *outdata,
			unsigned long int &outdataNcol, unsigned long int &outdataNrow,
			unsigned int narg, SEXP *argList) {
		if (indata) {
			outdata[0] = prod(indata, indataHeight*indataWidth);
		}
		outdataNcol = 1;
		outdataNrow = 1;
	}

	// Sum function + wrapper
	double sum(double *mydata, unsigned int size, bool dropNA) {
		double sumtotal = 0.;
		double zero = 0;
		//Rprintf("%f\n", mydata[0]);
		for (unsigned int i = 0; i < size; i++) {

			if (!ISNAN(mydata[i])) {
				sumtotal += mydata[i];
			} else if (!dropNA) {
				return(0/zero);
			}
		}
		return sumtotal;
	}

	void sumWrapper(double *indata, unsigned long int indataHeight,
			unsigned long int indataWidth,
			double *outdata,
			unsigned long int &outdataNcol, unsigned long int &outdataNrow,
			unsigned int narg, SEXP *argList) {
		if (indata) {
			bool dropNA = static_cast<bool> (INTEGER(argList[0])[0]);
			outdata[0] = sum(indata, indataHeight*indataWidth, dropNA);
		}
		outdataNcol = 1;
		outdataNrow = 1;
	}

	// Sum of powers function + wrapper
	double sumpower(double *mydata, unsigned int size, int power) {
		double sumpowertotal = 0.;
		for (unsigned int i = 0; i < size; i++) {
			sumpowertotal += pow(mydata[i], power);
		}
		return sumpowertotal;
	}
	void sumpowerWrapper(double *indata, unsigned long int indataHeight,
			unsigned long int indataWidth,
			double *outdata, unsigned long int &outdataNcol,
			unsigned long int &outdataNrow, unsigned int narg, SEXP *argList) {
		if (indata) {
			int power = INTEGER(argList[0])[0];
			outdata[0] = sumpower(indata, indataHeight*indataWidth, power);
		}
		outdataNcol = 1;
		outdataNrow = 1;
	}

//	void qtscore_globWrapper(double *indata, unsigned long int indataHeight,
//			unsigned long int indataWidth,
//			double *outdata, unsigned long int &outdataNcol,
//			unsigned long int &outdataNrow,	unsigned int narg, SEXP *argList) {
//		if(indata) {
//			// Fetching data from argList
//			double *pheno = REAL(argList[0]);
//			int    *Type  = INTEGER(argList[1]);
//			int    *Nids  = INTEGER(argList[2]);
//			int    *Nstra = INTEGER(argList[3]);
//			int    *stra  = INTEGER(argList[4]);
//			// Calling the qtscore_glob function
//			qtscore_glob(indata, pheno, Type, Nids, Nstra, stra, outdata);
//		}
//		outdataNcol = 10;
//		outdataNrow = 1;
//	}

//	void qtscore_glob(double *gt, double *pheno, int *Type, int *Nids, int *Nstra, int *stra, double *chi2)
//	{
//		int nsnps = 1;;
//		int nstra = (*Nstra);
//		int nids =  (*Nids);
//		int type = (*Type);
//		//int gt[nids];
//		int i, j, cstr, igt, i1=1;
//		int nbytes;
//		int dgt;
//		double Ttotg, mx, bb, totg[nstra], x2[nstra], sumx[nstra];
//		double Tsg0, Tsg1, Tsg2, sg0[nstra], sg1[nstra], sg2[nstra], xg0[nstra], xg1[nstra], xg2[nstra];
//		double u, v, u0, u1, u2, m0, m1, m2, v00, v02, v11, v12, v22, det;
//		mx = -999.99;
//		if ((nids % 4) == 0) nbytes = nids/4; else nbytes = ceil(1.*nids/4.);
//	//	char chgt[nbytes];
//
//		/*
//		 * The following loop has been disabled since iterator calls the qtscore_glob function
//		 * iterating over all snps.
//		 */
//		//for (igt=0;igt<nsnps;igt++) {
//		//	get_snps_many(gdata+nbytes*igt,Nids,&i1,gt);
//			for (j=0;j<nstra;j++) {
//				totg[j] = 0.;
//				x2[j] = 0.;
//				sumx[j] = 0.;
//				sg0[j] = 0.;
//				sg1[j] = 0.;
//				sg2[j] = 0.;
//				xg0[j] = 0.;
//				xg1[j] = 0.;
//				xg2[j] = 0.;
//			}
//
//			for (i=0;i<nids;i++) {
//			    if (gt[i] != 0) {
//			    	cstr = stra[i];
//			    	dgt = gt[i] - 1;
//			    	totg[cstr]+=1.0;
//			    	if (dgt==0) {
//			    		sg0[cstr]+=1.0;
//			    		xg0[cstr]+=pheno[i];
//			    	} else if (dgt==1) {
//			    		sg1[cstr]+=1.0;
//			    		xg1[cstr]+=pheno[i];
//			    	} else if (dgt==2) {
//			    		sg2[cstr]+=1.0;
//			    		xg2[cstr]+=pheno[i];
//			    	}
//			    	x2[cstr] += pheno[i]*pheno[i];
//			    	sumx[cstr] += pheno[i];
//			    }
//			}
//			Ttotg=Tsg0=Tsg1=Tsg2=0.;
//			for (j=0;j<nstra;j++) {
//				Ttotg += totg[j];
//				Tsg0 += sg0[j];
//				Tsg1 += sg1[j];
//				Tsg2 += sg2[j];
//			}
//			chi2[igt+6*nsnps]=Ttotg;
//			if (Ttotg == 0) {
//				chi2[igt] = -999.99;
//				chi2[igt+nsnps] = -999.99;
//				chi2[igt+2*nsnps] = -999.99;
//				chi2[igt+3*nsnps] = -999.99;
//				chi2[igt+4*nsnps] = -999.99;
//				chi2[igt+5*nsnps] = -999.99;
//				chi2[igt+7*nsnps] = -999.99;
//				chi2[igt+8*nsnps] = -999.99;
//				chi2[igt+9*nsnps] = -999.99;
//			} else {
//				u0 = u1 = u2 = m0 = m1 = m2 = v00 = v02 = v11 = v12 = v22 = 0.;
//				for (j=0;j<nstra;j++) if (totg[j]>0) {
//					mx = sumx[j]/totg[j];
//	//				if (type == 0)
//						bb = (x2[j]/totg[j])-mx*mx;
//	//				else
//	//					bb = mx*(1.-mx);
//					u0 += (xg0[j]-sg0[j]*mx);
//					m0 += xg0[j];
//					u1 += (xg1[j]-sg1[j]*mx);
//					m1 += xg1[j];
//					u2 += (xg2[j]-sg2[j]*mx);
//					m2 += xg2[j];
//					v00 += bb*(sg0[j]-sg0[j]*sg0[j]/totg[j]);
//					v11 += bb*(sg1[j]-sg1[j]*sg1[j]/totg[j]);
//					v12 += bb*(0.0-sg1[j]*sg2[j]/totg[j]);
//					v02 += bb*(0.0-sg0[j]*sg2[j]/totg[j]);
//					v22 += bb*(sg2[j]-sg2[j]*sg2[j]/totg[j]);
//				}
//				if (Tsg0>0) m0 = m0/Tsg0; else m0 = 1.e-16;
//				if (Tsg1>0) m1 = m1/Tsg1; else m1 = 1.e-16;
//				if (Tsg2>0) m2 = m2/Tsg2; else m2 = 1.e-16;
//				u = u1+2.*u2;
//				v = v11+4.*v12+4.*v22;
//				if (v<1.e-16) {
//				  chi2[igt]=-999.99;
//				  chi2[igt+3*nsnps]=-999.99;
//				} else {
//				  chi2[igt]=u*u/v;
//				  if (type) {
//					double p1 = mx+u/(Tsg1+4.*Tsg2-Ttotg*((Tsg1+2.*Tsg2)/Ttotg)*((Tsg1+2.*Tsg2)/Ttotg));
//				  	chi2[igt+3*nsnps]=(1.-mx)*p1/((1.-p1)*mx);
//				  } else {
//	//			  	chi2[igt+3*nsnps]=(Tsg0*(m0-mx)+Tsg1*(m1-mx)+Tsg2*(m2-mx))/Ttotg;
//				  	chi2[igt+3*nsnps]=u/(Tsg1+4.*Tsg2-Ttotg*((Tsg1+2.*Tsg2)/Ttotg)*((Tsg1+2.*Tsg2)/Ttotg));
//				  }
//				}
//				det = v11*v22 - v12*v12;
//	//			double rho2 = v12*v12/(v11*v22);
//	//			if (v00 <= 0. || v11<=0. || v22<=0. || rho2<1.e-16 || abs(det)<1.e-16) {
//					chi2[igt+nsnps] = -999.99;
//					chi2[igt+2*nsnps] = 1.e-16;
//				  	chi2[igt+4*nsnps] =-999.99;
//				  	chi2[igt+5*nsnps] = -999.99;
//					chi2[igt+7*nsnps] = -999.99;
//					chi2[igt+8*nsnps] = -999.99;
//					chi2[igt+9*nsnps] = -999.99;
//	//			} else {
//					if (v00>0.) {
//						chi2[igt+7*nsnps] = u0/sqrt(v00);
//						chi2[igt+nsnps] = u0*u0/v00;
//					}
//					if (v22>0.) {
//						chi2[igt+8*nsnps] = u2/sqrt(v22);
//						chi2[igt+nsnps] += u2*u2/v22;
//					}
//					if (v00*v22>0.) {
//						chi2[igt+9*nsnps] = v02/sqrt(v00*v22);
//						chi2[igt+nsnps] += -2.*u0*u2*v02/(v00*v22);
//						chi2[igt+nsnps] = chi2[igt+nsnps]/(1.-v02*v02/(v00*v22));
//					}
//	//				if (v11 > 0. && v22 > 0. && v12 > 0. && rho2<1.)
//	//				if (!(v00 <= 0. || v11<=0. || v22<=0. || rho2*rho2<1.e-16 || abs(det)<1.e-16))
//	//				HERE IS SOMETHING WRONG -- DO THE SAME AS IN QTSCORE CORRECTION!!!
//	//				if (!(v12 <= 0. || v11<=0. || v22<=0. || (rho2*rho2-1.)<1.e-16 || abs(det)<1.e-16))
//	//					chi2[igt+nsnps] = (u1*u1/v11 + u2*u2/v22 - 2.*rho2*u1*u2/v12)/(1.-rho2*rho2);
//	//				else
//	//					chi2[igt+nsnps] = chi2[igt];
//							//(u1*u1*v22+u2*u2*v11-2.0*u1*u2*v12)/det;
//					if (Tsg1>0) { if (type) {
//				  		chi2[igt+4*nsnps]=(1.-m0)*m1/((1.-m1)*m0);
//					} else {
//	//				 	chi2[igt+4*nsnps]=(u1/Tsg1);
//				  		chi2[igt+4*nsnps]=m1-m0;
//					}}
//					if (Tsg2>0) { if (type) {
//				  		chi2[igt+5*nsnps]=(1.-m0)*m2/((1.-m2)*m0);
//					} else {
//	//			  		chi2[igt+5*nsnps]=u2/Tsg2;
//				  		chi2[igt+5*nsnps]=m2-m0;
//					}}
//					if (Tsg1>0 && Tsg2>0)
//						chi2[igt+2*nsnps] = 2.;
//					else if (Tsg1>0 || Tsg2>0)
//						chi2[igt+2*nsnps] = 1.;
//	//			}
//			}
//	}


//void variance_homogeneity_test_C_wrapper(double *snp, unsigned long nids, double *outdata, 
//																				 unsigned long outdata_Ncol, unsigned long outdata_Nrow,
//	unsigned narg, double *argList)
//{
//}

//SEXP variance_homogeneity_test_C_wrapper_DEBUGING(SEXP trait_,
//																									SEXP design_matrix_,
//																									SEXP p_,
//																									SEXP analysis_type_,
//																									SEXP testname_,
//																									SEXP idnum_,
//																									SEXP snp_
//																									)
//{
//	
//
//	double *trait 	 			= REAL(trait_);
//	double *design_matrix = REAL(design_matrix_);
//	int    p  						= INTEGER_VALUE(p_); //covariate's number
//	int 	 analys_type    = INTEGER_VALUE(analysis_type_);
//	int    testname   	  = INTEGER_VALUE(testname_);
//	long unsigned int idnum = (long unsigned int)(INTEGER_VALUE(idnum_));
//	double *snp           = REAL(snp_);
//
//SEXP betas_;
//SEXP se_;
//SEXP chi2_;
//SEXP df_;
//SEXP residuals_;
//SEXP qty_;
//SEXP jpvt_;
//SEXP qraux_;
//SEXP work_;
//SEXP v_;
//SEXP x_for_ch2inv_;
//SEXP is_trait_na_;
//SEXP outdata_;
//
//
//PROTECT(betas_ = allocVector(REALSXP, p));
//PROTECT(se_ = allocVector(REALSXP, p));
//PROTECT(chi2_ = allocVector(REALSXP, 1));
//PROTECT(df_ = allocVector(INTSXP, 1));
//PROTECT(residuals_ = allocVector(REALSXP, idnum));
//PROTECT(qty_ = allocVector(REALSXP, idnum));
//PROTECT(jpvt_ = allocVector(INTSXP, p));
//PROTECT(qraux_ = allocVector(REALSXP, p));
//PROTECT(work_ = allocVector(REALSXP, p*2));
//PROTECT(v_ = allocVector(REALSXP, p*p));
//PROTECT(x_for_ch2inv_ = allocVector(REALSXP, p*p));
//PROTECT(is_trait_na_ = allocVector(INTSXP, idnum));
//PROTECT(outdata_ = allocVector(REALSXP, 2*(p) + 2));
//
//
//double *betas     		= REAL(betas_);
//double *se        		= REAL(se_);
//double *chi2      		= REAL(chi2_);
//int 	 *df        		= INTEGER(df_);
//double *residuals 		= REAL(residuals_);
///*auxiliary variables:*/
//double *qty       		= REAL(qty_);
//int *jpvt      		    = INTEGER(jpvt_);
//double *qraux     		= REAL(qraux_);
//double *work      		= REAL(work_);
//double *v         		= REAL(v_);
//double *x_for_ch2inv  = REAL(x_for_ch2inv_);
//int    *is_trait_na   = INTEGER(is_trait_na_);
////double *outdata       = REAL(outdata_);
//
//
//
//
//
//
//
//
//
////results_C <- .Call("variance_homogeneity_test_C_wrapper_DEBUGING",
////														 as.double(trait),
////														 as.double(data.matrix(design_matrix)),
////														 as.integer(p),
////														 as.integer(analysis_type),
////														 as.integer(testname),
////														 double(p),#betas 
////														 double(p),#se 
////														 double(1),#chi2
////														 integer(1),#df
////														 double(idnum),#residuals
////														 double(idnum),#qty 
////														 integer(p),#jpvt
////														 double(p),#qraux
////														 double(2*p),#work
////														 double(2*p),#v
////														 double(2*p),#x_for_ch2inv
////														 as.integer(idnum),
////														 integer(idnum),
////														 double(2*p + 2),
////														 as.double(snp)
////														 )
//
//
//
//
//
//
//
//	//Fill is_trait_na. 1 means missing vallue, 0 - non missing
//	//_________________________________________________________
//	for(unsigned id_counter=0 ; id_counter<idnum ; id_counter++)
//		{
//		if(ISNAN(trait[id_counter])) 
//			{
//			is_trait_na[id_counter]=1;
//			}
//		else
//			{
//			is_trait_na[id_counter]=0;
//			}
//		}
//	//_________________________________________________________
//
//
//
//	variance_homogeneity_test_C(snp, trait, design_matrix_geno_means, design_matrix, &p, &idnum, betas, se, chi2, df, residuals,
//				 											&analys_type, is_trait_na, &testname, /*auxiliary variables:*/ qty, jpvt, qraux, work, v, x_for_ch2inv);
//
//
//	REAL(outdata_)[0] = *chi2;
//	REAL(outdata_)[1] = double(*df);
//
//
//	for(int i=0 ; i<p ; i++)
//		{
//		REAL(outdata_)[i+2] = betas[i];
//		}
//	
//	for(int i=0 ; i<p ; i++)
//		{
//		REAL(outdata_)[i+2+p] = se[i];
//		}
//
//
////for(int i=0 ; i<(2+2*p) ; i++)
////	{
////	}			
//
//
//UNPROTECT(13);
//return outdata_;
//}
//



void variance_homogeneity_test_C_wrapper(double *indata, unsigned long int indataHeight,
																				 unsigned long int indataWidth,
																				 double *outdata, unsigned long int &outdataNcol,
																				 unsigned long int &outdataNrow,	unsigned int narg, SEXP *argList) 
{



//	static double *design_matrix_new = new double[indataHeight*(*p)];
//	static double *design_matrix_copy_new = new double[indataHeight*(*p)];

//	SEXP design_matrix_new_SEXP;
//	SEXP design_matrix_copy_new_SEXP;
//  PROTECT(design_matrix_new_SEXP = allocVector(REALSXP, indataHeight*(*p)));
//  PROTECT(design_matrix_copy_new_SEXP = allocVector(REALSXP, indataHeight*(*p)));
	
	
//	double *design_matrix_new = REAL(design_matrix_new_SEXP);
//	double *design_matrix_copy_new = REAL(design_matrix_copy_new_SEXP);
//	UNPROTECT(2);


//	for(int i=0 ; i<indataHeight*(*p) ; i++)
//		{
//		design_matrix_new[i] = design_matrix[i];
//		design_matrix_copy_new[i] = design_matrix_copy[i];
//		}


if(indata) 
	{
	double *trait_ 	 			           = REAL(argList[0]);
	double *design_matrix_            = REAL(argList[1]);
	double *design_matrix_copy_       = REAL(argList[2]);
	int    *p 											 = INTEGER(argList[3]); //covariate's number
	int 	 *analys_type              = INTEGER(argList[4]);
	int    *testname   	             = INTEGER(argList[5]);
	double *betas     		           = REAL(argList[6]);
	double *se        		           = REAL(argList[7]);
	double *chi2      		           = REAL(argList[8]);
	int 	 *df        	          	 = INTEGER(argList[9]);
	double *residuals 	 	           = REAL(argList[10]);
	/*auxiliary variables:*/
	double *qty       	 	           = REAL(argList[11]);
	int *jpvt      		               = INTEGER(argList[12]);
	double *qraux     		           = REAL(argList[13]);
	double *work      		           = REAL(argList[14]);
	double *v         		           = REAL(argList[15]);
	double *x_for_ch2inv             = REAL(argList[16]);
	int    *is_trait_na              = INTEGER(argList[17]);



	double *design_matrix = new double[indataHeight*(*p)];
	double *design_matrix_copy = new double[indataHeight*(*p)];
	double *trait = new double[indataHeight];	


	for(unsigned i=0 ; i<indataHeight*(*p) ; i++) design_matrix[i] = design_matrix_[i];
	for(unsigned i=0 ; i<indataHeight*(*p) ; i++) design_matrix_copy[i] = design_matrix_copy_[i];
	for(unsigned i=0 ; i<indataHeight ; i++) trait[i] = trait_[i];

	//indataHeight is nids;



	//Fill is_trait_na. 1 means missing vallue, 0 - non missing
	//_________________________________________________________
	for(unsigned id_counter=0 ; id_counter<indataHeight ; id_counter++)
		{
		is_trait_na[id_counter]=0;
		
		//NA among covariates
		for(int column_counter=0 ; column_counter<*p ; column_counter++)
			{
			if(ISNAN(design_matrix[id_counter + column_counter*indataHeight])) {is_trait_na[id_counter]=1;}
			}
	
		//Na among trait's values
		if(ISNAN(trait[id_counter])) 
			{
			is_trait_na[id_counter]=1;
			}
			
		//Na among SNP's values
		if(ISNAN(indata[id_counter])) 
			{
			is_trait_na[id_counter]=1;
			}
		

		}
	//_________________________________________________________


//	double *design_matrix_ = (double *) R_alloc(*p*indataWidth, sizeof(double));
//	double *design_matrix_ = new double[*p*indataWidth];
//  for(int i=0 ; i<*p*indataWidth ; i++) design_matrix_[i]=design_matrix[i];	

	


	variance_homogeneity_test_C(indata, trait, design_matrix_copy, design_matrix, p, &indataHeight, betas, se, chi2, df, residuals,
				 											analys_type, is_trait_na, testname, /*auxiliary variables:*/ qty, jpvt, qraux, work, v, x_for_ch2inv);

//	delete[] design_matrix_;

	outdata[0] = *chi2;
	outdata[1] = double(*df);

	int counter=2;
	for(int i=0 ; i<((*p)+1) ; i++)
		{
		outdata[counter] = betas[i];	
		counter++;
		outdata[counter] = se[i];	
		counter++;
		}

//	for(int i=0 ; i<2 + (*INTEGER(argList[3]))*2 ; i++)
//		{
//		}


	delete[] design_matrix;
	delete[] design_matrix_copy;
	delete[] trait;


	}
else
	{
	outdataNcol = 4 + (*INTEGER(argList[3]))*2;
	outdataNrow = 1;
	}
					
}

#ifdef __cplusplus
}
#endif
