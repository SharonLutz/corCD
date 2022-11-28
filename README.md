# corCD

The corCD R package checks whether the same direction of causality is returned when the regression models, from which the betas are pulled for each of the approaches, are adjusted for the covariates or unadjusted. It also checks the case returned and the absolute difference in K or the estimates between the adjusted and unadjusted models.


## Installation
```
install.packages("devtools") # devtools must be installed first
devtools::install_github("xue-hr/MRCD") # MRCD must be installed first
devtools::install_github("xue-hr/MRcML") # MRcML must be installed first
devtools::install_github("xue-hr/BiDirectCausal") # BiDirectCausal must be installed first

devtools::install_github("SharonLutz/corCD")
```

## Example
Under different scenarios given user input, one can test whether a method returns different cases depending on if the models to get the effect sizes adjusted for covariates or did not adjust for covariates. The below code runs an example.

```
library(corCD)
?corCD # for details on how to use this function

example1 <- corCD(n=70000, nSNPx = 20, MAFx = c(rep(0.2,  20)), nCovX = 2, meanCX = rep(0, 2), sdCX = c(rep(1, 2)),nSNPy = 20, MAFy = c(rep(0.2, 20)), nCovY = 2, meanCY = rep(0, 2), sdCY = c(rep(1,2)), deltaG = c(rep(0.2, 20)), deltaC =  c(rep(0.5,2),rep(1,2),rep(1.5,2),rep(2, 2),rep(2.5,2),rep(3,2)), betaG = c(rep(0.15, 20)), betaC=c(rep(0,2)), sdX=1,sdY=1, betaX=1, sig.level=0.05,  plot.pdf=F, table.sv=F, SEED = 1, nSims=50)
```

## Output
For this scenario, we have 5 matrices output. Four of the matrices return the number of simulations returning case 1 for each approaches, case 2 for each approach, case 3 for each approach, and the same case for adjusted and unadjusted models. One matrix returns the sum of the absolute difference of adjusted and unadjusted K or estimates.  

```
example1$mat1 # number of simulations returning case 1
example1$mat2 # number of simulations returning case 2
example1$mat3 # number of simulations returning case 3
example1$matSame # number of simulations returning the same case for adjusted and unadjusted models
example1$matDiff # sum of absolute difference between adjusted and unadjusted estimates
```

To get the average absolute difference in K or the estimates between the methods when using the unadjusted and adjusted models, as well as the proportion of simulations that return case 1, case 2, case 3, and the same case for the adjusted and unadjusted models, we can divide the matrices by the number of simulations run, in this example 50.

```
cbind(example1$mat1[,1], example1$mat1[,-1]/50) # proportion of simulations returning case 1
cbind(example1$mat2[,1], example1$mat2[,-1]/50) # proportion of simulations returning case 2
cbind(example1$mat3[,1], example1$mat3[,-1]/50) # proportion of simulations returning case 3
cbind(example1$matSame[,1], example1$matSame[,-1]/50) # proportion of simulations returning the same case for adjusted and unadjusted models
cbind(example1$matDiff[,1], example1$matDiff[,-1]/50) # average absolute difference between adjusted and unadjusted estimates
```    


