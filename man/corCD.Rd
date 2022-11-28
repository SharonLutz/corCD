\name{corCD}
\alias{corCD}
\title{corCD}
\description{The corCD function compares results from adjusted and unadjusted models when running functions from the following packages: "BiDirectCausal" and "MRCD". User input specifies the sample size (n), number of SNPs (nSNPx, nSNPy), minor allele frequencies (MAFx, MAFy), number of covariates (nCovX, nCovY), mean and standard deviation of the covariates (meanCX, meanCY, sdCX, sdCY), the effect sizes of the SNPs on X (deltaG), the effect sizes of the covariates on X (deltaC), the effect sizes of the SNPs on Y (betaG), the effect sizes of the covariates on Y (betaC), the standard deviation of X and Y (sdX, sdY), the effect size of X on Y (betaX), the seed (SEED), the significance level (sig.level), the number of simulations to run (nSims), whether to save the plots (plot.pdf), the name for the plots/tables (plot.name), and whether to save the tables (table.sv). The function is run for the user specified input, and each method is run when the models are (1) adjusted for covariates and (2) not adjusted for covariates.}
\usage{
corCD(n = 10000, nSNPx = 10, MAFx = c(rep(0.5, 10)), nCovX = 3, meanCX = rep(0, 3), sdCX = c(rep(1, 3)), nSNPy = 10, MAFy = c(rep(0.5, 10)), nCovY = 3, meanCY = rep(0, 3), sdCY = c(rep(1, 3)), deltaG = c(rep(0.2, 10)), deltaC = c(rep(0.7, 3), rep(0.9, 3), rep(1.1, 3)), betaG = c(rep(0.2, 10)), betaC = c(rep(0.5, 3)), sdX = 1, sdY = 1, betaX = 0.2, SEED = 1, sig.level = 0.05, nSims = 50,plot.pdf = T, plot.name = "plot", table.sv = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{is the sample size of the population that is being simulated.}
  \item{nSNPx}{is the number of SNPs for X.}
  \item{MAFx}{is the vector of minor allele frequencies for the SNPs for X.}
  \item{nCovX}{is the number of covariates for X.}
  \item{meanCX}{is the vector of the means for generating the SNPs for X.}
  \item{sdCX}{is the vector of standard deviations for generating the SNPs for X.}
  \item{nSNPy}{is the number of SNPs for Y.}
  \item{MAFy}{is the vector of minor allele frequencies for the SNPs for Y.}
  \item{nCovY}{is the number of covariates for Y.}
  \item{meanCY}{is the vector of the means for generating the SNPs for Y.}
  \item{sdCY}{is the vector of the standard deviations for generating the SNPs for Y.}
  \item{deltaG}{is the vector of effect sizes of the SNPs on X.}
  \item{deltaC}{is the vector of effect sizes of the covariates on X.}
  \item{betaG}{is the vector of effect sizes of the SNPs on Y.}
  \item{betaC}{is the vector of effect sizes of the covariates on Y.}
  \item{sdX}{is the standard deviation for generating X.}
  \item{sdY}{is the standard deviation for generating Y.}
  \item{betaX}{is the effect size of X on Y.}
  \item{SEED}{sets the seed used for the random generator.}
  \item{sig.level}{is the significance level.}
  \item{nSims}{is the number of simulations run for the function.}
  \item{plot.pdf}{if plot.pdf=T then 5 plots are saved to the working directory. If plot.pdf=F then no plots are created.}
  \item{plot.name}{specifies the name of the plots and tables that are created.}
  \item{table.sv}{if table.sv=T then the matrices output are also saved as .txt files in the working directory. If table.sv=F then the matrices are not saved as .txt files.}
}
\value{
    \item{mat1}{matrix of number of simulations that returned case 1 for each method.}
    \item{mat2}{matrix of number of simulations that returned case 2 for each method.}
    \item{mat3}{matrix of number of simulations that returned case 3 for each method.}
    \item{matDiff}{matrix of the sum of the absolute difference in K or estimate between the adjusted and unadjusted models.}
    \item{matSame}{matrix of the number of simulations that return the same case for the adjusted and unadjusted models.}
}
\references{
Xue H, Pan W. Inferring causal direction between two traits in the presence of horizontal pleiotropy with GWAS summary data. PLoS Genet. 2020 Nov 2;16(11):e1009105. doi: 10.1371/journal.pgen.1009105. PMID: 33137120; PMCID: PMC7660933.

Xue H, Pan W. Robust inference of bi-directional causal relationships in presence of correlated pleiotropy with GWAS summary data. PLoS Genet. 2022 May 16;18(5):e1010205. doi: 10.1371/journal.pgen.1010205. PMID: 35576237; PMCID: PMC9135345.
}
\author{
Sharon Lutz, Kirsten Voorhies
}

\section{Warning}{
library(MRCD), library(MRcML), and library(BiDirectCausal) are required to run this function.
}

\examples{test1 <- corCD(n=10000, nSNPx = 10, MAFx = c(rep(0.5,  10)), nCovX = 3, meanCX = rep(0, 3), sdCX = c(rep(1, 3)),
           nSNPy = 10, MAFy = c(rep(0.5, 10)), nCovY = 3, meanCY = rep(0, 3), sdCY = c(rep(1, 3)),
           deltaG = c(rep(0.2, 10)), deltaC =  c(rep(0.7,3), rep(0.9, 3), rep(1.1,3)), betaG = c(rep(0.2, 10)), betaC=c(rep(0.5,3)),
           sdX=1,sdY=1, betaX=0.2, SEED=1, sig.level=0.05, nSims=10, plot.pdf=FALSE, plot.name="plot", table.sv=FALSE)
           test1
           }

 \keyword{ function }
 \keyword{ correlation }

