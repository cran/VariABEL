\name{var_test_gwaa}
\alias{var_test_gwaa}
\title{Function performs comparison of genotypic variances.}
\description{The presence of a significant difference of genotypic variances
  points to possible presence of interaction between a tested SNP and a
  factor (or set of factors). This functions tests for that.
}

\usage{var_test_gwaa(formula, genodata, phenodata, genodata_info=NULL,
              testname="svlm", analysis_type="AAvsABvsBB")}

\arguments{
\item{formula}{Regression model used for analysis. In the first stage
  linear regression is run to exclude the main SNP effect. In this stage
  adjustment for covariates is performed.}

\item{genodata}{The genotypes data in GenABEL format or as DatABEL object.}

\item{phenodata}{The phenotypes data in format of data.frame object.}

\item{genodata_info}{The file with SNP information (name,
  position). Used if genotype data is a DatABEL object.}

\item{testname}{Name of the variance heterogeneity test to perform. "svlm"
  (for imputed genotype data), "Levene", and "Bartlett" test are supported. Note that "svlm" test only additive genetic model and the parameter analysis_type is ignored.}

\item{analysis_type}{Type of analysis to perform. \code{AAvsABvsBB} -
  each henotypic group is testes against other two, 
  \code{AAvsABandBB} - group AA tested against AB and BB,
  \code{ABvsAAandBB} - AB against AA and BB, \code{BBvsAAandAB} - BB
  against AA and AB. Only available for typed SNPs. Note that the input parameter analysis_type is ignored if testname="svlm"}
}

\details{
The function \code{var_test_gwaa} tests for difference in genotypic
variances. This difference points to the presence of possible interaction
between the tested SNP and some factor. In the case of the svlm test the
analysis consists of two stages: first the regular GWA is done where
regression analysis is performed with covariates specified in the input
parameter formula, in the second stage the regression analysis is
performed with using residuals from the first stage and a SNP as a
covariate.
}

\value{
The output is a \code{data.frame} object. The table contains the chi^2 of
the variance heterogeneity test (the name is \code{chisq}) the effects
and standard errors of all covariates included into the regression
model, main SNP effect (the names are \code{snp_eff} and
\code{snp_se}). In the case of the svlm test the columns
\code{snp_eff_dispertion} and \code{snp_se_dispertion} contain effect of
a SNP on the squared values of the trait.
}


%\references{ http://www.biomedcentral.com/1471-2156/11/92/abstract }
\author{Maksim Struchalin}
%\note{ ~~further notes~~

\examples{
if (require(GenABEL)) {
    data(srdta)
    result1 <- var_test_gwaa(bt ~ qt1 + qt2,
                             genodata=gtdata(srdta),
                             phenodata=phdata(srdta))

    # If there is covariates are needed:
    result2 <- var_test_gwaa("bt",
                             genodata=gtdata(srdta),
                             phenodata=phdata(srdta))
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
