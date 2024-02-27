# corCD

The corCD R package checks the direction of causality returned when the regression models, from which the betas are pulled for each of the approaches, are adjusted for the covariates or unadjusted. It also checks the difference in correlation of X and a SNP GX between the adjusted and unadjusted models as well as the difference in correlation of Y and a SNP GX between the adjusted and unadjusted models.

## Installation
```
install.packages("devtools") # devtools must be installed first
devtools::install_github("xue-hr/MRCD") # MRCD must be installed first
devtools::install_github("xue-hr/MRcML") # MRcML must be installed first
devtools::install_github("xue-hr/BiDirectCausal") # BiDirectCausal must be installed first
install.packages("RColorBrewer") # RColorBrewer must be installed first

devtools::install_github("SharonLutz/corCD")
```

## Example
Under different scenarios given user input, one can test which case a method returns depending on if the model effect sizes adjusted for covariates or did not adjust for covariates. The below code runs an example.

```
library(corCD)
?corCD # for details on how to use this function

example1 <- corCD(n = 1e+05, nSNPx = 10, MAFx = c(rep(0.5, 10)), nSNPy = 10, MAFy = c(rep(0.5, 10)),
nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 3)), deltaC = c(rep(0.005, 3)), betaC = c(rep(0.005, 3)),
betaGY = c(rep(0.2, 10)), betaGX = c(rep(0, 10)), deltaGX = c(rep(0.2, 10)), betaX = 0.2, sdX = 1,
sdY = 1, SEED = 1, sig.level = 0.05, nSims = 5, plot.pdf = TRUE, plot.name = "test", table.sv = TRUE,
approaches = c("All")) 
```

## Output
For this scenario, we have 3 matrices output. 1 matrix returns the number of simulations returning case 1 for each approach, case 2 for each approach, case 3 for each approach, for the adjusted and unadjusted models. One matrix returns the difference in correlation of cor(X, GX[,1]) when model effect sizes adjusted for covariates vs. did not adjust for covariates. One matrix returns the difference in correlation of cor(Y, GX[,1]) when model effect sizes adjusted for covariates vs. did not adjust for covariates

```
example1$matR # matrix of cases returned for methods run
example1$corGX # vector of difference in correlation of cor(X, GX[,1) when model effect sizes adjusted for covariates vs. did not adjust for covariates
example1$corGY # vector of difference in correlation of cor(Y, GX[,1) when model effect sizes adjusted for covariates vs. did not adjust for covariates
```

To get the proportion of simulations that return case 1, case 2, case 3 for each method methods when using the unadjusted and adjusted models, we can divide the matrix example1$matR by the number of simulations run, in this example 10.

```
cbind(example1$matR/50) # proportion of simulations returning the different cases for the different methods

cbind(example1$matSame[,1], example1$matSame[,-1]/50) # proportion of simulations returning the same case for adjusted and unadjusted models
cbind(example1$matDiff[,1], example1$matDiff[,-1]/50) # average absolute difference between adjusted and unadjusted estimates
```    


