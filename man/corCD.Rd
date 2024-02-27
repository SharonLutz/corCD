\name{corCD}
\alias{corCD}
\title{corCD}
\description{
The function runs simulations and outputs (1) the proportion of simulations returning the different cases for the simple linear regression and multiple linear regression for the approaches, (2) the difference in correlation between X and the first SNP in the matrix GX for the multiple linear regression and simple linear regression, (3) the difference in correlation between Y and the first SNP in the matrix GX for the multiple linear regression and simple linear regression.
}
\usage{
corCD(n = 1e+05, nSNPx = 10, MAFx = c(rep(0.5, 10)), nSNPy = 10, MAFy = c(rep(0.5, 10)), nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 3)), deltaC = c(rep(0.005, 3)), betaC = c(rep(0.005, 3)), betaGY = c(rep(0.2, 10)), betaGX = c(rep(0, 10)), deltaGX = c(rep(0.2, 10)), betaX = 0.2, sdX = 1, sdY = 1, SEED = 1, sig.level = 0.05, nSims = 5, plot.pdf = FALSE, plot.name = "plot", table.sv = FALSE, approaches = c("All"))
}

\arguments{
  \item{n}{
is the sample size of the population that is being simulated.
}
  \item{nSNPx}{
is the number of SNPs for X.
}
  \item{MAFx}{
is the vector of minor allele frequencies for the SNPs for X.
}
  \item{nSNPy}{
is the number of SNPs for Y.
}
  \item{MAFy}{
is the vector of minor allele frequencies for the SNPs for Y.
}
  \item{nCov}{
is the number of covariates.
}
  \item{meanC}{
is the vector of the means for generating the covariates.
}
  \item{sdC}{
is the vector of the standard deviations for generating the covariates.
}
  \item{deltaC}{
is the vector of effect sizes of the covariates on X.
}
  \item{betaC}{
is the vector of effect sizes of the covariates on Y.
}
  \item{betaGY}{
is the vector of effect sizes of the SNPs GY on Y.
}
  \item{betaGX}{
is the vector of effect sizes of the SNPs GX on Y.
}
  \item{deltaGX}{
is the vector of effect sizes of the SNPs GX on X.
}
  \item{betaX}{
is the effect size of X on Y.
}
  \item{sdX}{
is the standard deviation for generating X.
 }
  \item{sdY}{
is the standard deviation for generating Y.
}
  \item{SEED}{
sets the seed used for the random generator.
}
  \item{sig.level}{
is the significance level.
}
  \item{nSims}{
is the number of simulations run for the function.
}
  \item{plot.pdf}{
if plot.pdf=TRUE, plots of the results will be saved as PDFs to the working directory.
}
  \item{plot.name}{
specifies the name to use for saving the plots and tables.
}
  \item{table.sv}{
if table.sv=TRUE, a table of the results will be saved as a .txt file to the working directory.
}
  \item{approaches}{
specifies which approaches to run. Default is "All", which will run CDcML, MRcML, CD-Ratio, CD-GLS, and CD-Egger. 
One or more of these methods can be selected by providing a list of one or more of the following: "CDcML", "MRcML", "CDRatio", "CDEgger", "CDGLS".
}
}

\value{
\item{matR}{matrix of cases returned for each of the approaches run.}
\item{corGX}{vector of the difference in correlation between the multiple linear regression and simple linear regression when X is the outcome and the first SNP in GX is the SNP.}
\item{corGY}{vector of the difference in correlation between the multiple linear regression and simple linear regression when Y is the outcome and the first SNP in GX is the SNP.}
}
\references{
Xue H, Pan W. Inferring causal direction between two traits in the presence of horizontal pleiotropy with GWAS summary data. PLoS Genet. 2020 Nov 2;16(11):e1009105. doi: 10.1371/journal.pgen.1009105. PMID: 33137120; PMCID: PMC7660933.

Xue H, Pan W. Robust inference of bi-directional causal relationships in presence of correlated pleiotropy with GWAS summary data. PLoS Genet. 2022 May 16;18(5):e1010205. doi: 10.1371/journal.pgen.1010205. PMID: 35576237; PMCID: PMC9135345.
}
\author{
Sharon Lutz, Kirsten Voorhies
}
\section{Warning }{
library(MRCD), library(MRcML), library(BiDirectCausal), and library(RColorBrewer) are needed to run this function.
}
\examples{
testC <- corCD(n = 1e+05, nSNPx = 10, MAFx = c(rep(0.5, 10)), nSNPy = 10, 
    MAFy = c(rep(0.5, 10)), nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 
        3)), deltaC = c(rep(0.005, 3)), betaC = c(rep(0.005, 
        3)), betaGY = c(rep(0.2, 10)), betaGX = c(rep(0, 10)), 
    deltaGX = c(rep(0.2, 10)), betaX = 0.2, sdX = 1, sdY = 1, 
    SEED = 1, sig.level = 0.05, nSims = 5, plot.pdf = TRUE, plot.name = "test", 
    table.sv = TRUE, approaches = c("All")) 
}
\keyword{ function }
\keyword{ correlation }

