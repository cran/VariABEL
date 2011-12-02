\name{var_test_gwaa}
\alias{var_test_gwaa}
\title{Function performs compariosn of genotypic variances.}
\description{
		Presence of significant difference of genotypic variances points to possible 
						Significant variance difference points to possible presence of interaction between a tested SNP and a factor (or set of factors).
}
\usage{
		var_test_gwaa(formula, genodata, phenodata, genodata_info=NULL, testname="svlm", analysis_type="AAvsABvsBB")
}
\arguments{\item{formula}{Regression model used for analysis. In the first stage linear regression is run to exclude main snp effect. In this stage adjustment for covariates is performed.}
		\item{genodata}{The genotypes data in format of genabel or databel object.}
		\item{phenodata}{The phenotypes data in format of data.frame object.}
		\item{genodata_info}{The file with snp information (name, position). Used if genodata is databel object.}
		\item{testname}{Name of variance heterogeneity test to perform. svlm (for imputed genotype data), levene, and bartlett test are supported.}
		\item{analysis_type}{Type of analsysis to perform. AAvsABvsBB - additive model where B allele additivly increase risk, AAvsABandBB - group AA tested agains AB and BB, ABvsAAandBB - AB against AA and BB, BBvsAAandAB - BB against AA and AB. Only available for typed snps.}
}
\details{

		The function var_test_gwaa tests for difference in genotypic variances. This difference points to presence of possible interaction between the tested SNP 
				and some factor. In the case svlm test the analysis consists of two stage: firstly the regular GWA id done where regression analysis is performed with 
				covariates specified in the input parameter formula, in the second stage the regression analysis is performed with using residuals from the first stage and 
				a sno as a covariate.
}

\value{
			The ouput is a data.frame object. The table contains the chisq of variance heterogeneity test (the name is chisq) the effects and standart errors of all covariates included into regression model, 
						main snp effect (the names are snp_eff and snp_se). In the case of svlm test the columns snp_eff_dispertion and snp_se_dispertion contain
									effect of a snp on squared vallues of a trait.

										}
  

%\references{ http://www.biomedcentral.com/1471-2156/11/92/abstract }
\author{Maksim Struchalin}
%\note{ ~~further notes~~ 
\examples{
	if (require(GenABEL))
 		{
		data(srdta)
		result1 <- var_test_gwaa(bt~qt1+qt2, genodata=gtdata(srdta), phenodata=phdata(srdta))
														
#if there is no covariates needed:
		result2 <- var_test_gwaa("bt", genodata=gtdata(srdta), phenodata=phdata(srdta))
		}

	}

\references{

Struchalin et al., Variance heterogeneity analysis for detection of potentially 
interacting genetic loci: method and its limitations. 
BMC Genetics 2010, 11:92, doi:10.1186/1471-2156-11-92
								
Struchalin et al., An R package VariABEL for genome-wide searching of
potentially interacting loci by testing genotypic variance heterogeneity. Submitted. 
										
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}


