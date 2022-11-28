corCD <-
function(n=10000, nSNPx = 10, MAFx = c(rep(0.5,  10)), nCovX = 3, meanCX = rep(0, 3), sdCX = c(rep(1, 3)),
           nSNPy = 10, MAFy = c(rep(0.5, 10)), nCovY = 3, meanCY = rep(0, 3), sdCY = c(rep(1, 3)),
           deltaG = c(rep(0.2, 10)), deltaC =  c(rep(0.7,3), rep(0.9, 3), rep(1.1,3)), betaG = c(rep(0.2, 10)), betaC=c(rep(0.5,3)),
           sdX=1,sdY=1, betaX=0.2, SEED=1, sig.level=0.05, nSims=50, plot.pdf=T, plot.name="plot", table.sv=F){
    
    
    deltaC<-matrix(deltaC, ncol=nCovX, byrow=T)
    
    if(length(MAFx)!=nSNPx){stop("length(MAFx) must equal nSNPx.")}
    if(length(MAFy)!=nSNPy){stop("length(MAFy) must equal nSNPy")}
    if(length(deltaG)!=nSNPx){stop("length(deltaG) must equal nSNPx")}
    if(length(betaG)!=nSNPy){stop("length(betaG) must equal nSNPy")}
    
    if(length(meanCX)!=nCovX){stop("length(meanCX) must equal nCovX")}
    if(length(meanCY)!=nCovY){stop("length(meanCY) must equal nCovY")}
    if(length(sdCX)!=nCovX){stop("length(sdCX) must equal nCovX")}
    if(length(sdCY)!=nCovY){stop("length(sdCY) must equal nCovY")}
    if(ncol(deltaC)!=nCovX){stop("ncol(deltaC) must equal nCovX")}
    if(length(betaC)!=nCovY){stop("length(betaC) must equal nCovY")}
    
    # load libraries
    library(MRCD)
    library(MRcML)
    library(BiDirectCausal)
    ###################################
    # create matrices, assign column names
    mat1 <- matrix(0, ncol=31, nrow=nrow(deltaC))
    mat2 <- matrix(0, ncol=31, nrow=nrow(deltaC))
    mat3 <- matrix(0, ncol=31, nrow=nrow(deltaC))
    matDiff <- matrix(0, ncol=31, nrow=nrow(deltaC))
    matSame <- matrix(0, ncol=16, nrow = nrow(deltaC))
    
    colnames.use1 <- c("deltaC","CDcML.S.DP.XY", "CDcML.S.noDP.XY", "CDcML.noS.DP.XY", "CDcML.noS.noDP.XY",
                       "CDcML.S.DP.YX", "CDcML.S.noDP.YX", "CDcML.noS.DP.YX", "CDcML.noS.noDP.YX",
                       "MRcML.S.DP.XY", "MRcML.S.noDP.XY", "MRcML.noS.DP.XY", "MRcML.noS.noDP.XY",
                       "MRcML.S.DP.YX", "MRcML.S.noDP.YX", "MRcML.noS.DP.YX", "MRcML.noS.noDP.YX",
                       "biCDRatioXY.S", "biCDRatioXY.noS", "biCDEggerXY.S", "biCDEggerXY.noS",
                       "biCDRatioYX.S", "biCDRatioYX.noS", "biCDEggerYX.S", "biCDEggerYX.noS",
                       "CDRatioXY", "CDEggerXY", "CDGLSXY",
                       "CDRatioYX", "CDEggerYX", "CDGLSYX")
    
    colnames.use2 <- c("deltaC", "CDcML.S.DP.u", "CDcML.S.noDP.u", "CDcML.noS.DP.u", "CDcML.noS.noDP.u",
                       "CDcML.S.DP.a", "CDcML.S.noDP.a", "CDcML.noS.DP.a", "CDcML.noS.noDP.a",
                       "MRcML.S.DP.u", "MRcML.S.noDP.u", "MRcML.noS.DP.u", "MRcML.noS.noDP.u",
                       "MRcML.S.DP.a", "MRcML.S.noDP.a", "MRcML.noS.DP.a", "MRcML.noS.noDP.a",
                       "biCDRatio.S.u", "biCDRatio.noS.u", "biCDEgger.S.u", "biCDEgger.noS.u",
                       "biCDRatio.S.a", "biCDRatio.noS.a", "biCDEgger.S.a", "biCDEgger.noS.a",
                       "CDRatio.u", "CDEgger.u", "CDGLS.u",
                       "CDRatio.a", "CDEgger.a", "CDGLS.a")
    
    colnames.use3 <- c("deltaC", "CDcML.S.DP", "CDcML.S.noDP", "CDcML.noS.DP", "CDcML.noS.noDP",
                       "MRcML.S.DP", "MRcML.S.noDP", "MRcML.noS.DP", "MRcML.noS.noDP",
                       "biCDRatio.S", "biCDRatio.noS", "biCDEgger.S", "biCDEgger.noS",
                       "CDRatio", "CDEgger", "CDGLS")
    
    colnames(mat1) <- colnames.use2
    colnames(mat2) <- colnames.use2
    colnames(mat3) <- colnames.use2
    colnames(matDiff) <- colnames.use1
    colnames(matSame) <- colnames.use3
    
    mat1[,"deltaC"] <- deltaC[,1]
    mat2[,"deltaC"] <- deltaC[,1]
    mat3[,"deltaC"] <- deltaC[,1]
    matDiff[,"deltaC"] <- deltaC[,1]
    matSame[,"deltaC"] <- deltaC[,1]
    
    
    for(j in 1:nrow(deltaC)){  
      
      for(i in 1:nSims){
        printCut<-1
        if(floor(i/printCut)==ceiling(i/printCut)){print(paste(i,"of",nSims, "simulations and", j, "of",nrow(deltaC),"deltaC"))}
        set.seed(SEED+i-1)
        
        ###################################
        # matrix for deltaG, deltaC, betaG
        deltaG<-matrix(deltaG, ncol=1, nrow=nSNPx)
        betaC<-matrix(betaC, ncol=1, nrow=nCovY)
        betaG<- matrix(betaG,ncol=1,nrow=nSNPy)
        
        # create SNP matrix for X
        matGX<-matrix(0,ncol=nSNPx,nrow=n)
        for(ss in 1:nSNPx){
          matGX[,ss]<-rbinom(n,2, MAFx[ss])
        }
        
        # create SNP matrix for Y
        matGY<-matrix(0,ncol=nSNPy,nrow=n)
        for(ss in 1:nSNPy){
          matGY[,ss]<-rbinom(n,2, MAFy[ss])
        }
        
        # create Covariate matrix
        matCX<-matrix(0, ncol=nCovX,nrow=n)
        for(cc in 1:nCovX){
          matCX[,cc]<-rnorm(n, meanCX[cc], sdCX[cc])
        }
        
        # create Covariate matrix
        matCY<-matrix(0, ncol=nCovY,nrow=n)
        for(cc in 1:nCovY){
          matCY[,cc]<-rnorm(n, meanCY[cc], sdCY[cc])
        }
        x<-rnorm(n, (matGX%*%deltaG+ matCX%*%deltaC[j,]), sdX) # trait X
        y<- rnorm(n, (matGY%*%betaG+ matCY%*%betaC+betaX*x), sdY) # trait X
        
        ###################################
        # GX --> X --> Y
        bxU <- NULL
        sxU <- NULL
        byU <- NULL
        syU <- NULL
        for(ll in 1:ncol(matGX)){
          bxU <- c(bxU, summary(lm(x~matGX[,ll]))$coef[2,1])
          sxU <- c(sxU, summary(lm(x~matGX[,ll]))$coef[2,2])
          
          byU <- c(byU, summary(lm(y~matGX[,ll]))$coef[2,1])
          syU <- c(syU, summary(lm(y~matGX[,ll]))$coef[2,2])
        }
        
        for(ll in 1:ncol(matGY)){
          bxU <- c(bxU, summary(lm(x~matGY[,ll]))$coef[2,1])
          sxU <- c(sxU, summary(lm(x~matGY[,ll]))$coef[2,2])
          
          byU <- c(byU, summary(lm(y~matGY[,ll]))$coef[2,1])
          syU <- c(syU, summary(lm(y~matGY[,ll]))$coef[2,2])
        }
        
        bxA <- NULL
        sxA <- NULL
        byA <- NULL
        syA <- NULL
        for(ll in 1:ncol(matGX)){
          bxA <- c(bxA, summary(lm(x~matGX[,ll]+matCX))$coef[2,1])
          sxA <- c(sxA, summary(lm(x~matGX[,ll]+matCX))$coef[2,2])
          
          byA <- c(byA, summary(lm(y~matGX[,ll]+matCY))$coef[2,1])
          syA <- c(syA, summary(lm(y~matGX[,ll]+matCY))$coef[2,2])
        }
        
        for(ll in 1:ncol(matGY)){
          bxA <- c(bxA, summary(lm(x~matGY[,ll]+matCX))$coef[2,1])
          sxA <- c(sxA, summary(lm(x~matGY[,ll]+matCX))$coef[2,2])
          
          byA <- c(byA, summary(lm(y~matGY[,ll]+matCY))$coef[2,1])
          syA <- c(syA, summary(lm(y~matGY[,ll]+matCY))$coef[2,2])
        }
        
        ###################################
        # CDcML
        res.BDCDcML.u <- BiDirCDcML(b_X = bxU,
                                    b_Y = byU,
                                    se_X = sxU,
                                    se_Y = syU,
                                    n_X = n,
                                    n_Y = n,
                                    sig.cutoff = 0.05/(nSNPx+nSNPy),
                                    num_pert = 100)
        
        # CDcML
        res.BDCDcML.a <- BiDirCDcML(b_X = bxA,
                                    b_Y = byA,
                                    se_X = sxA,
                                    se_Y = syA,
                                    n_X = n,
                                    n_Y = n,
                                    sig.cutoff = 0.05/(nSNPx+nSNPy),
                                    num_pert = 100)
        ############################################
        # Unadjusted - CDcCML
        ############################################
        
        u.S.DP<-NA
        u.S.noDP<-NA
        u.noS.DP<-NA
        u.noS.noDP<-NA
        # CIs - .S.DP
        u.ci.low.XY.S.DP <- res.BDCDcML.u$XtoY.est.S.DP - (qnorm(0.975)*res.BDCDcML.u$XtoY.se.S.DP)
        u.ci.upp.XY.S.DP <- res.BDCDcML.u$XtoY.est.S.DP + (qnorm(0.975)*res.BDCDcML.u$XtoY.se.S.DP)
        
        u.ci.low.YX.S.DP <- res.BDCDcML.u$YtoX.est.S.DP - (qnorm(0.975)*res.BDCDcML.u$YtoX.se.S.DP)
        u.ci.upp.YX.S.DP <- res.BDCDcML.u$YtoX.est.S.DP + (qnorm(0.975)*res.BDCDcML.u$YtoX.se.S.DP)
        
        # decision rules with CIs  --- .S.DP
        if((u.ci.low.XY.S.DP>0 & u.ci.upp.XY.S.DP<1)| (u.ci.low.XY.S.DP> (-1) & u.ci.upp.XY.S.DP<0)){
          if((u.ci.low.YX.S.DP>0 & u.ci.upp.YX.S.DP<1) | (u.ci.low.YX.S.DP> (-1) & u.ci.upp.YX.S.DP<0)){
            mat3[j,"CDcML.S.DP.u"] <- mat3[j,"CDcML.S.DP.u"]+1
            u.S.DP <- 3
          }else{
            mat1[j,"CDcML.S.DP.u"] <- mat1[j,"CDcML.S.DP.u"]+1
            u.S.DP <- 1
          }
        }else if((u.ci.low.YX.S.DP>0 & u.ci.upp.YX.S.DP<1) & (u.ci.low.YX.S.DP> (-1) & u.ci.upp.YX.S.DP<0)){
          if((u.ci.low.XY.S.DP>0 & u.ci.upp.XY.S.DP<1)| (u.ci.low.XY.S.DP> (-1) & u.ci.upp.XY.S.DP<0)){
            mat3[j,"CDcML.S.DP.u"] <- mat3[j,"CDcML.S.DP.u"]+1
            u.S.DP <- 3
          }else{mat2[j,"CDcML.S.DP.u"] <- mat2[j,"CDcML.S.DP.u"]+1
          u.S.DP <- 2
          }
        }else{mat3[j,"CDcML.S.DP.u"] <- mat3[j,"CDcML.S.DP.u"]+1
        u.S.DP <- 3
        }
        
        # CIs - NoS
        u.ci.low.XY.NoS <- res.BDCDcML.u$XtoY.est.NoS - (qnorm(0.975)*res.BDCDcML.u$XtoY.se.NoS)
        u.ci.upp.XY.NoS <- res.BDCDcML.u$XtoY.est.NoS + (qnorm(0.975)*res.BDCDcML.u$XtoY.se.NoS)
        
        u.ci.low.YX.NoS <- res.BDCDcML.u$YtoX.est.NoS - (qnorm(0.975)*res.BDCDcML.u$YtoX.se.NoS)
        u.ci.upp.YX.NoS <- res.BDCDcML.u$YtoX.est.NoS + (qnorm(0.975)*res.BDCDcML.u$YtoX.se.NoS)        
        # decision rules with CIs  --- NoS
        if((u.ci.low.XY.NoS>0 & u.ci.upp.XY.NoS<1)| (u.ci.low.XY.NoS>(-1) & u.ci.upp.XY.NoS<0)){
          if((u.ci.low.YX.NoS>0 & u.ci.upp.YX.NoS<1) | (u.ci.low.YX.NoS> (-1) & u.ci.upp.YX.NoS<0)){
            mat3[j,"CDcML.noS.noDP.u"] <- mat3[j,"CDcML.noS.noDP.u"]+1
            u.noS.noDP <- 3
          }else{
            mat1[j,"CDcML.noS.noDP.u"] <- mat1[j,"CDcML.noS.noDP.u"]+1
            u.noS.noDP <- 1
          }
        }else if((u.ci.low.YX.NoS>0 & u.ci.upp.YX.NoS<1) & (u.ci.low.YX.NoS>(-1) & u.ci.upp.YX.NoS<0)){
          if((u.ci.low.XY.NoS>0 & u.ci.upp.XY.NoS<1)| (u.ci.low.XY.NoS>(-1) & u.ci.upp.XY.NoS<0)){
            mat3[j,"CDcML.noS.noDP.u"] <- mat3[j,"CDcML.noS.noDP.u"]+1
            u.noS.noDP <- 3
          }else{
            mat2[j,"CDcML.noS.noDP.u"] <- mat2[j,"CDcML.noS.noDP.u"]+1
            u.noS.noDP <- 2
          }
        }else{
          mat3[j,"CDcML.noS.noDP.u"] <- mat3[j,"CDcML.noS.noDP.u"]+1
          u.noS.noDP <- 3
        }
        
        # CIs - S
        u.ci.low.XY.S <- res.BDCDcML.u$XtoY.est.S - (qnorm(0.975)*res.BDCDcML.u$XtoY.se.S)
        u.ci.upp.XY.S <- res.BDCDcML.u$XtoY.est.S + (qnorm(0.975)*res.BDCDcML.u$XtoY.se.S)
        
        u.ci.low.YX.S <- res.BDCDcML.u$YtoX.est.S - (qnorm(0.975)*res.BDCDcML.u$YtoX.se.S)
        u.ci.upp.YX.S <- res.BDCDcML.u$YtoX.est.S + (qnorm(0.975)*res.BDCDcML.u$YtoX.se.S)  
        # decision rules with CIs  --- S
        if((u.ci.low.XY.S>0 & u.ci.upp.XY.S<1)| (u.ci.low.XY.S>-1 & u.ci.upp.XY.S<0)){
          if((u.ci.low.YX.S>0 & u.ci.upp.YX.S<1) | (u.ci.low.YX.S>-1 & u.ci.upp.YX.S<0)){
            mat3[j,"CDcML.S.noDP.u"] <- mat3[j,"CDcML.S.noDP.u"]+1
            u.S.noDP <- 3
          }else{
            mat1[j,"CDcML.S.noDP.u"] <- mat1[j,"CDcML.S.noDP.u"]+1
            u.S.noDP <- 1
          }
        }else if((u.ci.low.YX.S>0 & u.ci.upp.YX.S<1) & (u.ci.low.YX.S>-1 & u.ci.upp.YX.S<0)){
          if((u.ci.low.XY.S>0 & u.ci.upp.XY.S<1)| (u.ci.low.XY.S>-1 & u.ci.upp.XY.S<0)){
            mat3[j,"CDcML.S.noDP.u"] <- mat3[j,"CDcML.S.noDP.u"]+1
            u.S.noDP <- 3
          }else{
            mat2[j,"CDcML.S.noDP.u"] <- mat2[j,"CDcML.S.noDP.u"]+1
            u.S.noDP <- 2
          }
        }else{
          mat3[j,"CDcML.S.noDP.u"] <- mat3[j,"CDcML.S.noDP.u"]+1
          u.S.noDP <- 3
        }
        
        
        # CIs - NoS.DP
        u.ci.low.XY.NoS.DP <- res.BDCDcML.u$XtoY.est.NoS.DP - (qnorm(0.975)*res.BDCDcML.u$XtoY.se.NoS.DP)
        u.ci.upp.XY.NoS.DP <- res.BDCDcML.u$XtoY.est.NoS.DP + (qnorm(0.975)*res.BDCDcML.u$XtoY.se.NoS.DP)
        
        u.ci.low.YX.NoS.DP <- res.BDCDcML.u$YtoX.est.NoS.DP - (qnorm(0.975)*res.BDCDcML.u$YtoX.se.NoS.DP)
        u.ci.upp.YX.NoS.DP <- res.BDCDcML.u$YtoX.est.NoS.DP + (qnorm(0.975)*res.BDCDcML.u$YtoX.se.NoS.DP)        
        # decision rules with CIs  --- NoS
        if((u.ci.low.XY.NoS.DP>0 & u.ci.upp.XY.NoS.DP<1)| (u.ci.low.XY.NoS.DP>-1 & u.ci.upp.XY.NoS.DP<0)){
          if((u.ci.low.YX.NoS.DP>0 & u.ci.upp.YX.NoS.DP<1) | (u.ci.low.YX.NoS.DP>-1 & u.ci.upp.YX.NoS.DP<0)){
            mat3[j,"CDcML.noS.DP.u"] <- mat3[j,"CDcML.noS.DP.u"]+1
            u.noS.DP <- 3
          }else{
            mat1[j,"CDcML.noS.DP.u"] <- mat1[j,"CDcML.noS.DP.u"]+1
            u.noS.DP <- 1
          }
        }else if((u.ci.low.YX.NoS.DP>0 & u.ci.upp.YX.NoS.DP<1) & (u.ci.low.YX.NoS.DP>-1 & u.ci.upp.YX.NoS.DP<0)){
          if((u.ci.low.XY.NoS.DP>0 & u.ci.upp.XY.NoS.DP<1)| (u.ci.low.XY.NoS.DP>-1 & u.ci.upp.XY.NoS.DP<0)){
            mat3[j,"CDcML.noS.DP.u"] <- mat3[j,"CDcML.noS.DP.u"]+1
            u.noS.DP <- 3
          }else{
            mat2[j,"CDcML.noS.DP.u"] <- mat2[j,"CDcML.noS.DP.u"]+1
            u.noS.DP <- 2
          }
        }else{
          mat3[j,"CDcML.noS.DP.u"] <- mat3[j,"CDcML.noS.DP.u"]+1
          u.noS.DP <- 3
        }
        
        ############################################
        # Adjusted - CDcML
        ############################################
        a.S.DP<-NA
        a.S.noDP<-NA
        a.noS.DP<-NA
        a.noS.noDP<-NA
        # CIs - .S.DP
        a.ci.low.XY.S.DP <- res.BDCDcML.a$XtoY.est.S.DP - (qnorm(0.975)*res.BDCDcML.a$XtoY.se.S.DP)
        a.ci.upp.XY.S.DP <- res.BDCDcML.a$XtoY.est.S.DP + (qnorm(0.975)*res.BDCDcML.a$XtoY.se.S.DP)
        
        a.ci.low.YX.S.DP <- res.BDCDcML.a$YtoX.est.S.DP - (qnorm(0.975)*res.BDCDcML.a$YtoX.se.S.DP)
        a.ci.upp.YX.S.DP <- res.BDCDcML.a$YtoX.est.S.DP + (qnorm(0.975)*res.BDCDcML.a$YtoX.se.S.DP)
        
        # decision rules with CIs  --- .S.DP
        if((a.ci.low.XY.S.DP>0 & a.ci.upp.XY.S.DP<1)| (a.ci.low.XY.S.DP>-1 & a.ci.upp.XY.S.DP<0)){
          if((a.ci.low.YX.S.DP>0 & a.ci.upp.YX.S.DP<1) | (a.ci.low.YX.S.DP>-1 & a.ci.upp.YX.S.DP<0)){
            mat3[j,"CDcML.S.DP.a"] <- mat3[j,"CDcML.S.DP.a"]+1
            a.S.DP<- 3
          }else{
            mat1[j,"CDcML.S.DP.a"] <- mat1[j,"CDcML.S.DP.a"]+1
            a.S.DP<- 1
          }
        }else if((a.ci.low.YX.S.DP>0 & a.ci.upp.YX.S.DP<1) & (a.ci.low.YX.S.DP>-1 & a.ci.upp.YX.S.DP<0)){
          if((a.ci.low.XY.S.DP>0 & a.ci.upp.XY.S.DP<1)| (a.ci.low.XY.S.DP>-1 & a.ci.upp.XY.S.DP<0)){
            mat3[j,"CDcML.S.DP.a"] <- mat3[j,"CDcML.S.DP.a"]+1
            a.S.DP<- 3
          }else{
            mat2[j,"CDcML.S.DP.a"] <- mat2[j,"CDcML.S.DP.a"]+1
            a.S.DP<- 2
          }
        }else{
          mat3[j,"CDcML.S.DP.a"] <- mat3[j,"CDcML.S.DP.a"]+1
          a.S.DP<- 3
        }
        
        # CIs - NoS
        a.ci.low.XY.NoS <- res.BDCDcML.a$XtoY.est.NoS - (qnorm(0.975)*res.BDCDcML.a$XtoY.se.NoS)
        a.ci.upp.XY.NoS <- res.BDCDcML.a$XtoY.est.NoS + (qnorm(0.975)*res.BDCDcML.a$XtoY.se.NoS)
        
        a.ci.low.YX.NoS <- res.BDCDcML.a$YtoX.est.NoS - (qnorm(0.975)*res.BDCDcML.a$YtoX.se.NoS)
        a.ci.upp.YX.NoS <- res.BDCDcML.a$YtoX.est.NoS + (qnorm(0.975)*res.BDCDcML.a$YtoX.se.NoS)        
        # decision rules with CIs  --- NoS
        if((a.ci.low.XY.NoS>0 & a.ci.upp.XY.NoS<1)| (a.ci.low.XY.NoS>-1 & a.ci.upp.XY.NoS<0)){
          if((a.ci.low.YX.NoS>0 & a.ci.upp.YX.NoS<1) | (a.ci.low.YX.NoS>-1 & a.ci.upp.YX.NoS<0)){
            mat3[j,"CDcML.noS.noDP.a"] <- mat3[j,"CDcML.noS.noDP.a"]+1
            a.noS.noDP<- 3
          }else{
            mat1[j,"CDcML.noS.noDP.a"] <- mat1[j,"CDcML.noS.noDP.a"]+1
            a.noS.noDP<- 1
          }
        }else if((a.ci.low.YX.NoS>0 & a.ci.upp.YX.NoS<1) & (a.ci.low.YX.NoS>-1 & a.ci.upp.YX.NoS<0)){
          if((a.ci.low.XY.NoS>0 & a.ci.upp.XY.NoS<1)| (a.ci.low.XY.NoS>-1 & a.ci.upp.XY.NoS<0)){
            mat3[j,"CDcML.noS.noDP.a"] <- mat3[j,"CDcML.noS.noDP.a"]+1
            a.noS.noDP<- 3
          }else{
            mat2[j,"CDcML.noS.noDP.a"] <- mat2[j,"CDcML.noS.noDP.a"]+1
            a.noS.noDP<- 2
          }
        }else{
          mat3[j,"CDcML.noS.noDP.a"] <- mat3[j,"CDcML.noS.noDP.a"]+1
          a.noS.noDP<- 3
        }
        
        # CIs - Screening, no DP
        a.ci.low.XY.S <- res.BDCDcML.a$XtoY.est.S - (qnorm(0.975)*res.BDCDcML.a$XtoY.se.S)
        a.ci.upp.XY.S <- res.BDCDcML.a$XtoY.est.S + (qnorm(0.975)*res.BDCDcML.a$XtoY.se.S)
        
        a.ci.low.YX.S <- res.BDCDcML.a$YtoX.est.S - (qnorm(0.975)*res.BDCDcML.a$YtoX.se.S)
        a.ci.upp.YX.S <- res.BDCDcML.a$YtoX.est.S + (qnorm(0.975)*res.BDCDcML.a$YtoX.se.S)  
        # decision rules with CIs  --- S
        if((a.ci.low.XY.S>0 & a.ci.upp.XY.S<1)| (a.ci.low.XY.S>-1 & a.ci.upp.XY.S<0)){
          if((a.ci.low.YX.S>0 & a.ci.upp.YX.S<1) | (a.ci.low.YX.S>-1 & a.ci.upp.YX.S<0)){
            mat3[j,"CDcML.S.noDP.a"] <- mat3[j,"CDcML.S.noDP.a"]+1
            a.S.noDP<- 3
          }else{
            mat1[j,"CDcML.S.noDP.a"] <- mat1[j,"CDcML.S.noDP.a"]+1
            a.S.noDP<- 1
          }
        }else if((a.ci.low.YX.S>0 & a.ci.upp.YX.S<1) & (a.ci.low.YX.S>-1 & a.ci.upp.YX.S<0)){
          if((a.ci.low.XY.S>0 & a.ci.upp.XY.S<1)| (a.ci.low.XY.S>-1 & a.ci.upp.XY.S<0)){
            mat3[j,"CDcML.S.noDP.a"] <- mat3[j,"CDcML.S.noDP.a"]+1
            a.S.noDP<- 3
          }else{
            mat2[j,"CDcML.S.noDP.a"] <- mat2[j,"CDcML.S.noDP.a"]+1
            a.S.noDP<- 2
          }
        }else{
          mat3[j,"CDcML.S.noDP.a"] <- mat3[j,"CDcML.S.noDP.a"]+1
          a.S.noDP<- 3
        }
        
        
        # CIs - No screening; yes DP
        a.ci.low.XY.NoS.DP <- res.BDCDcML.a$XtoY.est.NoS.DP - (qnorm(0.975)*res.BDCDcML.a$XtoY.se.NoS.DP)
        a.ci.upp.XY.NoS.DP <- res.BDCDcML.a$XtoY.est.NoS.DP + (qnorm(0.975)*res.BDCDcML.a$XtoY.se.NoS.DP)
        
        a.ci.low.YX.NoS.DP <- res.BDCDcML.a$YtoX.est.NoS.DP - (qnorm(0.975)*res.BDCDcML.a$YtoX.se.NoS.DP)
        a.ci.upp.YX.NoS.DP <- res.BDCDcML.a$YtoX.est.NoS.DP + (qnorm(0.975)*res.BDCDcML.a$YtoX.se.NoS.DP)        
        # decision rules with CIs  --- NoS
        if((a.ci.low.XY.NoS.DP>0 & a.ci.upp.XY.NoS.DP<1)| (a.ci.low.XY.NoS.DP>-1 & a.ci.upp.XY.NoS.DP<0)){
          if((a.ci.low.YX.NoS.DP>0 & a.ci.upp.YX.NoS.DP<1) | (a.ci.low.YX.NoS.DP>-1 & a.ci.upp.YX.NoS.DP<0)){
            mat3[j,"CDcML.noS.DP.a"] <- mat3[j,"CDcML.noS.DP.a"]+1
            a.noS.DP <- 3
          }else{
            mat1[j,"CDcML.noS.DP.a"] <- mat1[j,"CDcML.noS.DP.a"]+1
            a.noS.DP <- 1
          }
        }else if((a.ci.low.YX.NoS.DP>0 & a.ci.upp.YX.NoS.DP<1) & (a.ci.low.YX.NoS.DP>-1 & a.ci.upp.YX.NoS.DP<0)){
          if((a.ci.low.XY.NoS.DP>0 & a.ci.upp.XY.NoS.DP<1)| (a.ci.low.XY.NoS.DP>-1 & a.ci.upp.XY.NoS.DP<0)){
            mat3[j,"CDcML.noS.DP.a"] <- mat3[j,"CDcML.noS.DP.a"]+1
            a.noS.DP <- 3
          }else{
            mat2[j,"CDcML.noS.DP.a"] <- mat2[j,"CDcML.noS.DP.a"]+1
            a.noS.DP <- 2
          }
        }else{
          mat3[j,"CDcML.noS.DP.a"] <- mat3[j,"CDcML.noS.DP.a"]+1
          a.noS.DP <- 3
        }
        
        # fill in matDiff for CD s and no s approaches
        matDiff[j,"CDcML.S.DP.XY"] <- matDiff[j,"CDcML.S.DP.XY"] + abs(res.BDCDcML.a$XtoY.est.S.DP - res.BDCDcML.u$XtoY.est.S.DP)
        matDiff[j,"CDcML.S.DP.YX"] <- matDiff[j,"CDcML.S.DP.YX"]  + abs(res.BDCDcML.a$YtoX.est.S.DP - res.BDCDcML.u$YtoX.est.S.DP)
        
        matDiff[j,"CDcML.S.noDP.XY"] <- matDiff[j,"CDcML.S.noDP.XY"] + abs(res.BDCDcML.a$XtoY.est.S - res.BDCDcML.u$XtoY.est.S)
        matDiff[j,"CDcML.S.noDP.YX"] <- matDiff[j,"CDcML.S.noDP.YX"] + abs(res.BDCDcML.a$YtoX.est.S - res.BDCDcML.u$YtoX.est.S)
        
        matDiff[j,"CDcML.noS.DP.XY"] <- matDiff[j,"CDcML.noS.DP.XY"] + abs(res.BDCDcML.a$XtoY.est.NoS.DP - res.BDCDcML.u$XtoY.est.NoS.DP)
        matDiff[j,"CDcML.noS.DP.YX"] <- matDiff[j,"CDcML.noS.DP.YX"] + abs(res.BDCDcML.a$YtoX.est.NoS.DP - res.BDCDcML.u$YtoX.est.NoS.DP)
        
        matDiff[j,"CDcML.noS.noDP.XY"] <- matDiff[j,"CDcML.noS.noDP.XY"] + abs(res.BDCDcML.a$XtoY.est.NoS - res.BDCDcML.u$XtoY.est.NoS)
        matDiff[j,"CDcML.noS.noDP.YX"] <- matDiff[j,"CDcML.noS.noDP.YX"] + abs(res.BDCDcML.a$YtoX.est.NoS - res.BDCDcML.u$YtoX.est.NoS)
        
        
        # Fill in matSame for CDcML
        
        if(u.S.DP==a.S.DP){matSame[j,"CDcML.S.DP"] <- matSame[j,"CDcML.S.DP"] + 1}
        
        if(u.S.noDP==a.S.noDP){matSame[j,"CDcML.S.noDP"] <- matSame[j,"CDcML.S.noDP"] + 1}
        
        if(u.noS.DP==a.noS.DP){matSame[j,"CDcML.noS.DP"] <- matSame[j,"CDcML.noS.DP"] + 1}
        
        if(u.noS.noDP==a.noS.noDP){matSame[j,"CDcML.noS.noDP"] <- matSame[j,"CDcML.noS.noDP"] + 1}
        
        ##########################################
        # MRcML
        ##########################################
        res.BDMRcML.u <-BiDirMRcML(b_X = bxU,
                                   b_Y = byU,
                                   se_X = sxU,
                                   se_Y = syU,
                                   n_X = n,
                                   n_Y = n,
                                   sig.cutoff = 0.05/(nSNPx+nSNPy),
                                   num_pert = 100)
        res.BDMRcML.a <-BiDirMRcML(b_X = bxA,
                                   b_Y = byA,
                                   se_X = sxA,
                                   se_Y = syA,
                                   n_X = n,
                                   n_Y = n,
                                   sig.cutoff = 0.05/(nSNPx+nSNPy),
                                   num_pert = 100)
        
        ############################################
        # Unadjusted - MRcCML
        ############################################
        u.S.DP<-NA
        u.S.noDP<-NA
        u.noS.DP<-NA
        u.noS.noDP<-NA
        
        u.p.XY.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.u$XtoY.est.S.DP / res.BDMRcML.u$XtoY.se.S.DP))*2
        u.p.YX.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.u$YtoX.est.S.DP / res.BDMRcML.u$YtoX.se.S.DP))*2
        # decision rules
        if(u.p.XY.BDMRcML.S.DP <= 0.05 & u.p.YX.BDMRcML.S.DP > 0.05){
          mat1[j,"MRcML.S.DP.u"] <- mat1[j,"MRcML.S.DP.u"]+1
          u.S.DP<-1
        }#return case 1
        if(u.p.XY.BDMRcML.S.DP > 0.05 & u.p.YX.BDMRcML.S.DP <= 0.05){
          mat2[j,"MRcML.S.DP.u"] <- mat2[j,"MRcML.S.DP.u"]+1
          u.S.DP<-2
        }# return case 2
        if((u.p.XY.BDMRcML.S.DP >0.05 & u.p.YX.BDMRcML.S.DP > 0.05) |(u.p.XY.BDMRcML.S.DP <=0.05 & u.p.YX.BDMRcML.S.DP <= 0.05)){
          mat3[j,"MRcML.S.DP.u"] <- mat3[j,"MRcML.S.DP.u"]+1
          u.S.DP<-3
        }# return case 3
        
        
        # P-value bi-directional MRcML methods: .S
        u.p.XY.BDMRcML.S<- pnorm(-abs(res.BDMRcML.u$XtoY.est.S / res.BDMRcML.u$XtoY.se.S))*2
        u.p.YX.BDMRcML.S<- pnorm(-abs(res.BDMRcML.u$YtoX.est.S / res.BDMRcML.u$YtoX.se.S))*2
        # decision rules
        if(u.p.XY.BDMRcML.S <= 0.05 & u.p.YX.BDMRcML.S > 0.05){
          mat1[j,"MRcML.S.noDP.u"] <- mat1[j,"MRcML.S.noDP.u"]+1
          u.S.noDP<-1
        }#return case 1
        if(u.p.XY.BDMRcML.S > 0.05 & u.p.YX.BDMRcML.S <= 0.05){
          mat2[j,"MRcML.S.noDP.u"] <- mat2[j,"MRcML.S.noDP.u"]+1
          u.S.noDP<-2
        }# return case 2
        if((u.p.XY.BDMRcML.S >0.05 & u.p.YX.BDMRcML.S > 0.05) |(u.p.XY.BDMRcML.S <=0.05 & u.p.YX.BDMRcML.S <= 0.05)){
          mat3[j,"MRcML.S.noDP.u"] <- mat3[j,"MRcML.S.noDP.u"]+1
          u.S.noDP<-3
        }# return case 3
        
        
        # P-value bi-directional MRcML methods: .NoS
        u.p.XY.BDMRcML.NoS<- pnorm(-abs(res.BDMRcML.u$XtoY.est.NoS / res.BDMRcML.u$XtoY.se.NoS))*2
        u.p.YX.BDMRcML.NoS<- pnorm(-abs(res.BDMRcML.u$YtoX.est.NoS / res.BDMRcML.u$YtoX.se.NoS))*2
        # decision rules
        if(u.p.XY.BDMRcML.NoS <= 0.05 & u.p.YX.BDMRcML.NoS > 0.05){
          mat1[j,"MRcML.noS.noDP.u"] <- mat1[j,"MRcML.noS.noDP.u"]+1
          u.noS.noDP<-1
        }#return case 1
        if(u.p.XY.BDMRcML.NoS > 0.05 & u.p.YX.BDMRcML.NoS <= 0.05){
          mat2[j,"MRcML.noS.noDP.u"] <- mat2[j,"MRcML.noS.noDP.u"]+1
          u.noS.noDP<-2
        }# return case 2
        if((u.p.XY.BDMRcML.NoS >0.05 & u.p.YX.BDMRcML.NoS > 0.05) |(u.p.XY.BDMRcML.NoS <=0.05 & u.p.YX.BDMRcML.NoS <= 0.05)){
          mat3[j,"MRcML.noS.noDP.u"] <- mat3[j,"MRcML.noS.noDP.u"]+1
          u.noS.noDP<-3
        }# return case 3
        
        # P-value bi-directional MRcML methods: .NoS.DP
        u.p.XY.BDMRcML.NoS.DP<- pnorm(-abs(res.BDMRcML.u$XtoY.est.NoS.DP / res.BDMRcML.u$XtoY.se.NoS.DP))*2
        u.p.YX.BDMRcML.NoS.DP<- pnorm(-abs(res.BDMRcML.u$YtoX.est.NoS.DP / res.BDMRcML.u$YtoX.se.NoS.DP))*2
        # decision rules
        if(u.p.XY.BDMRcML.NoS.DP <= 0.05 & u.p.YX.BDMRcML.NoS.DP > 0.05){
          mat1[j,"MRcML.noS.DP.u"] <- mat1[j,"MRcML.noS.DP.u"]+1
          u.noS.DP<-1
        }#return case 1
        if(u.p.XY.BDMRcML.NoS.DP > 0.05 & u.p.YX.BDMRcML.NoS.DP <= 0.05){
          mat2[j,"MRcML.noS.DP.u"] <- mat2[j,"MRcML.noS.DP.u"]+1
          u.noS.DP<-2
        }# return case 2
        if((u.p.XY.BDMRcML.NoS.DP >0.05 & u.p.YX.BDMRcML.NoS.DP > 0.05) |(u.p.XY.BDMRcML.NoS.DP <=0.05 & u.p.YX.BDMRcML.NoS.DP <= 0.05)){
          mat3[j,"MRcML.noS.DP.u"] <- mat3[j,"MRcML.noS.DP.u"]+1
          u.noS.DP<-3
        }# return case 3
        
        ############################################
        # Adjusted - MRcCML
        ############################################
        a.S.DP<-NA
        a.S.noDP<-NA
        a.noS.DP<-NA
        a.noS.noDP<-NA
        
        a.p.XY.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.a$XtoY.est.S.DP / res.BDMRcML.a$XtoY.se.S.DP))*2
        a.p.YX.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.a$YtoX.est.S.DP / res.BDMRcML.a$YtoX.se.S.DP))*2
        # decision rules
        if(a.p.XY.BDMRcML.S.DP <= 0.05 & a.p.YX.BDMRcML.S.DP > 0.05){
          mat1[j,"MRcML.S.DP.a"] <- mat1[j,"MRcML.S.DP.a"]+1
          a.S.DP<-1
        }#return case 1
        if(a.p.XY.BDMRcML.S.DP > 0.05 & a.p.YX.BDMRcML.S.DP <= 0.05){
          mat2[j,"MRcML.S.DP.a"] <- mat2[j,"MRcML.S.DP.a"]+1
          a.S.DP<-2
        }# return case 2
        if((a.p.XY.BDMRcML.S.DP >0.05 & a.p.YX.BDMRcML.S.DP > 0.05) |(a.p.XY.BDMRcML.S.DP <=0.05 & a.p.YX.BDMRcML.S.DP <= 0.05)){
          mat3[j,"MRcML.S.DP.a"] <- mat3[j,"MRcML.S.DP.a"]+1
          a.S.DP<-3
        }# return case 3
        
        # P-value bi-directional MRcML methods: .S
        a.p.XY.BDMRcML.S<- pnorm(-abs(res.BDMRcML.a$XtoY.est.S / res.BDMRcML.a$XtoY.se.S))*2
        a.p.YX.BDMRcML.S<- pnorm(-abs(res.BDMRcML.a$YtoX.est.S / res.BDMRcML.a$YtoX.se.S))*2
        # decision rules
        if(a.p.XY.BDMRcML.S <= 0.05 & a.p.YX.BDMRcML.S > 0.05){
          mat1[j,"MRcML.S.noDP.a"] <- mat1[j,"MRcML.S.noDP.a"]+1
          a.S.noDP<-1
        }#return case 1
        if(a.p.XY.BDMRcML.S > 0.05 & a.p.YX.BDMRcML.S <= 0.05){
          mat2[j,"MRcML.S.noDP.a"] <- mat2[j,"MRcML.S.noDP.a"]+1
          a.S.noDP<-2
        }# return case 2
        if((a.p.XY.BDMRcML.S >0.05 & a.p.YX.BDMRcML.S > 0.05) |(a.p.XY.BDMRcML.S <=0.05 & a.p.YX.BDMRcML.S <= 0.05)){
          mat3[j,"MRcML.S.noDP.a"] <- mat3[j,"MRcML.S.noDP.a"]+1
          a.S.noDP<-3
        }# return case 3
        
        # P-value bi-directional MRcML methods: .NoS
        a.p.XY.BDMRcML.NoS<- pnorm(-abs(res.BDMRcML.a$XtoY.est.NoS / res.BDMRcML.a$XtoY.se.NoS))*2
        a.p.YX.BDMRcML.NoS<- pnorm(-abs(res.BDMRcML.a$YtoX.est.NoS / res.BDMRcML.a$YtoX.se.NoS))*2
        # decision rules
        if(a.p.XY.BDMRcML.NoS <= 0.05 & a.p.YX.BDMRcML.NoS > 0.05){
          mat1[j,"MRcML.noS.noDP.a"] <- mat1[j,"MRcML.noS.noDP.a"]+1
          a.noS.noDP<-1
        }#return case 1
        if(a.p.XY.BDMRcML.NoS > 0.05 & a.p.YX.BDMRcML.NoS <= 0.05){
          mat2[j,"MRcML.noS.noDP.a"] <- mat2[j,"MRcML.noS.noDP.a"]+1
          a.noS.noDP<-2
        }# return case 2
        if((a.p.XY.BDMRcML.NoS >0.05 & a.p.YX.BDMRcML.NoS > 0.05) |(a.p.XY.BDMRcML.NoS <=0.05 & a.p.YX.BDMRcML.NoS <= 0.05)){
          mat3[j,"MRcML.noS.noDP.a"] <- mat3[j,"MRcML.noS.noDP.a"]+1
          a.noS.noDP<-3
        }# return case 3
        
        # P-value bi-directional MRcML methods: .NoS.DP
        a.p.XY.BDMRcML.NoS.DP<- pnorm(-abs(res.BDMRcML.a$XtoY.est.NoS.DP / res.BDMRcML.a$XtoY.se.NoS.DP))*2
        a.p.YX.BDMRcML.NoS.DP<- pnorm(-abs(res.BDMRcML.a$YtoX.est.NoS.DP / res.BDMRcML.a$YtoX.se.NoS.DP))*2
        # decision rules
        if(a.p.XY.BDMRcML.NoS.DP <= 0.05 & a.p.YX.BDMRcML.NoS.DP > 0.05){
          mat1[j,"MRcML.noS.DP.a"] <- mat1[j,"MRcML.noS.DP.a"]+1
          a.noS.DP<-1
        }#return case 1
        if(a.p.XY.BDMRcML.NoS.DP > 0.05 & a.p.YX.BDMRcML.NoS.DP <= 0.05){
          mat2[j,"MRcML.noS.DP.a"] <- mat2[j,"MRcML.noS.DP.a"]+1
          a.noS.DP<-2
        }# return case 2
        if((a.p.XY.BDMRcML.NoS.DP >0.05 & a.p.YX.BDMRcML.NoS.DP > 0.05) |(a.p.XY.BDMRcML.NoS.DP <=0.05 & a.p.YX.BDMRcML.NoS.DP <= 0.05)){
          mat3[j,"MRcML.noS.DP.a"] <- mat3[j,"MRcML.noS.DP.a"]+1
          a.noS.DP<-3
        }# return case 3
        
        ##################################################################
        # fill in matDiff for MRcML approaches
        ##################################################################
        matDiff[j,"MRcML.S.DP.XY"] <- matDiff[j,"MRcML.S.DP.XY"] + abs(res.BDMRcML.a$XtoY.est.S.DP - res.BDMRcML.u$XtoY.est.S.DP)
        matDiff[j,"MRcML.S.DP.YX"] <- matDiff[j,"MRcML.S.DP.YX"] + abs(res.BDMRcML.a$YtoX.est.S.DP - res.BDMRcML.u$YtoX.est.S.DP)
        
        matDiff[j,"MRcML.S.noDP.XY"] <- matDiff[j,"MRcML.S.noDP.XY"] + abs(res.BDMRcML.a$XtoY.est.S - res.BDMRcML.u$XtoY.est.S)
        matDiff[j,"MRcML.S.noDP.YX"] <- matDiff[j,"MRcML.S.noDP.YX"] + abs(res.BDMRcML.a$YtoX.est.S - res.BDMRcML.u$YtoX.est.S)
        
        matDiff[j,"MRcML.noS.DP.XY"] <- matDiff[j,"MRcML.noS.DP.XY"] + abs(res.BDMRcML.a$XtoY.est.NoS.DP - res.BDMRcML.u$XtoY.est.NoS.DP)
        matDiff[j,"MRcML.noS.DP.YX"] <- matDiff[j,"MRcML.noS.DP.YX"] + abs(res.BDMRcML.a$YtoX.est.NoS.DP - res.BDMRcML.u$YtoX.est.NoS.DP)
        
        matDiff[j,"MRcML.noS.noDP.XY"] <- matDiff[j,"MRcML.noS.noDP.XY"] + abs(res.BDMRcML.a$XtoY.est.NoS - res.BDMRcML.u$XtoY.est.NoS)
        matDiff[j,"MRcML.noS.noDP.YX"] <- matDiff[j,"MRcML.noS.noDP.YX"] + abs(res.BDMRcML.a$YtoX.est.NoS - res.BDMRcML.u$YtoX.est.NoS)
        
        
        ##################################################################
        # Fill in matSame for MRcML
        ##################################################################
        if(u.S.DP==a.S.DP){matSame[j,"MRcML.S.DP"] <- matSame[j,"MRcML.S.DP"] + 1}
        
        if(u.S.noDP==a.S.noDP){matSame[j,"MRcML.S.noDP"] <- matSame[j,"MRcML.S.noDP"] + 1}
        
        if(u.noS.DP==a.noS.DP){matSame[j,"MRcML.noS.DP"] <- matSame[j,"MRcML.noS.DP"] + 1}
        
        if(u.noS.noDP==a.noS.noDP){matSame[j,"MRcML.noS.noDP"] <- matSame[j,"MRcML.noS.noDP"] + 1}
        
        
        ################################################
        # Run BiDirCDMethod
        ################################################
        # CD-Ratio - s &  CD-Egger - s - unadjusted
        res.biCD.U <- BiDirCDMethod(b_X = bxU,
                                    b_Y = byU,
                                    se_X = sxU,
                                    se_Y = syU,
                                    n_X = n,
                                    n_Y = n,
                                    sig.cutoff = 0.05/(nSNPx+nSNPy))
        
        
        # CD-Ratio - s &  CD-Egger - s - adjusted
        res.biCD.A <- BiDirCDMethod(b_X = bxA,
                                    b_Y = byA,
                                    se_X = sxA,
                                    se_Y = syA,
                                    n_X = n,
                                    n_Y = n,
                                    sig.cutoff = 0.05/(nSNPx+nSNPy))
        
        ##########################################
        # Results for unadjusted with screening
        ##########################################
        
        uR.S <- NA
        uR.noS <- NA
        uE.S <- NA
        uE.noS <- NA
        
        u.K.ratio.XY.S <- res.biCD.U$CDRatio.XtoY.est.S
        u.seK.ratio.XY.S <- res.biCD.U$CDRatio.XtoY.se.S
        u.K.ratio.YX.S <- res.biCD.U$CDRatio.YtoX.est.S
        u.seK.ratio.YX.S <- res.biCD.U$CDRatio.YtoX.se.S
        
        u.K.egger.XY.S <- res.biCD.U$CDEgger.XtoY.est.S
        u.seK.egger.XY.S <- res.biCD.U$CDEgger.XtoY.se.S
        u.K.egger.YX.S <- res.biCD.U$CDEgger.YtoX.est.S
        u.seK.egger.YX.S <- res.biCD.U$CDEgger.YtoX.se.S
        
        u.ci.low.ratioXY.S <- u.K.ratio.XY.S - (qnorm(0.975)*u.seK.ratio.XY.S)
        u.ci.upp.ratioXY.S <- u.K.ratio.XY.S + (qnorm(0.975)*u.seK.ratio.XY.S)
        u.ci.low.ratioYX.S <- u.K.ratio.YX.S - (qnorm(0.975)*u.seK.ratio.YX.S)
        u.ci.upp.ratioYX.S <- u.K.ratio.YX.S + (qnorm(0.975)*u.seK.ratio.YX.S)
        
        u.ci.low.eggerXY.S <- u.K.egger.XY.S - (qnorm(0.975)*u.seK.egger.XY.S)
        u.ci.upp.eggerXY.S <- u.K.egger.XY.S + (qnorm(0.975)*u.seK.egger.XY.S)
        u.ci.low.eggerYX.S <- u.K.egger.YX.S - (qnorm(0.975)*u.seK.egger.YX.S)
        u.ci.upp.eggerYX.S <- u.K.egger.YX.S + (qnorm(0.975)*u.seK.egger.YX.S)
        
        #  DECISION RULES CD Ratio - .S unadjusted
        if((u.ci.low.ratioXY.S>0 & u.ci.upp.ratioXY.S<1)| (u.ci.low.ratioXY.S>(-1) & u.ci.upp.ratioXY.S<0)){
          if((u.ci.low.ratioYX.S>0 & u.ci.upp.ratioYX.S<1) | (u.ci.low.ratioYX.S>(-1) & u.ci.upp.ratioYX.S<0)){
            mat3[j,"biCDRatio.S.u"] <- mat3[j,"biCDRatio.S.u"]+1
            uR.S <- 3
          }else{
            mat1[j,"biCDRatio.S.u"] <- mat1[j,"biCDRatio.S.u"]+1
            uR.S <- 1
          }
        }else if((u.ci.low.ratioYX.S>0 & u.ci.upp.ratioYX.S<1) & (u.ci.low.ratioYX.S>(-1) & u.ci.upp.ratioYX.S<0)){
          if((u.ci.low.ratioXY.S>0 & u.ci.upp.ratioXY.S<1)| (u.ci.low.ratioXY.S>(-1) & u.ci.upp.ratioXY.S<0)){
            mat3[j,"biCDRatio.S.u"] <- mat3[j,"biCDRatio.S.u"]+1
            uR.S <- 3
          }else{
            mat2[j,"biCDRatio.S.u"] <- mat2[j,"biCDRatio.S.u"]+1
            uR.S <- 2
          }
        }else{
          mat3[j,"biCDRatio.S.u"] <- mat3[j,"biCDRatio.S.u"]+1
          uR.S <- 3
        }
        
        #  DECISION RULES CD Egger - .S unadjusted
        if((u.ci.low.eggerXY.S>0 & u.ci.upp.eggerXY.S<1)| (u.ci.low.eggerXY.S>(-1) & u.ci.upp.eggerXY.S<0)){
          if((u.ci.low.eggerYX.S>0 & u.ci.upp.eggerYX.S<1) | (u.ci.low.eggerYX.S>(-1) & u.ci.upp.eggerYX.S<0)){
            mat3[j,"biCDEgger.S.u"] <- mat3[j,"biCDEgger.S.u"]+1
            uE.S <- 3
          }else{mat1[j,"biCDEgger.S.u"] <- mat1[j,"biCDEgger.S.u"]+1
          uE.S <- 1
          }
        }else if((u.ci.low.eggerYX.S>0 & u.ci.upp.eggerYX.S<1) & (u.ci.low.eggerYX.S>(-1) & u.ci.upp.eggerYX.S<0)){
          if((u.ci.low.eggerXY.S>0 & u.ci.upp.eggerXY.S<1)| (u.ci.low.eggerXY.S>(-1) & u.ci.upp.eggerXY.S<0)){
            mat3[j,"biCDEgger.S.u"] <- mat3[j,"biCDEgger.S.u"]+1
            uE.S <- 3
          }else{
            mat2[j,"biCDEgger.S.u"] <- mat2[j,"biCDEgger.S.u"]+1
            uE.S <- 2
          }
        }else{
          mat3[j,"biCDEgger.S.u"] <- mat3[j,"biCDEgger.S.u"]+1
          uE.S <- 3
        }
        
        # save K & SE for no screening - unadjusted
        u.K.ratio.XY.NoS <- res.biCD.U$CDRatio.XtoY.est.NoS
        u.seK.ratio.XY.NoS <- res.biCD.U$CDRatio.XtoY.se.NoS
        u.K.ratio.YX.NoS <- res.biCD.U$CDRatio.YtoX.est.NoS
        u.seK.ratio.YX.NoS <- res.biCD.U$CDRatio.YtoX.se.NoS
        
        u.K.egger.XY.NoS <- res.biCD.U$CDEgger.XtoY.est.NoS
        u.seK.egger.XY.NoS <- res.biCD.U$CDEgger.XtoY.se.NoS
        u.K.egger.YX.NoS <- res.biCD.U$CDEgger.YtoX.est.NoS
        u.seK.egger.YX.NoS <- res.biCD.U$CDEgger.YtoX.se.NoS
        
        u.ci.low.ratioXY.NoS <- u.K.ratio.XY.NoS - (qnorm(0.975)*u.seK.ratio.XY.NoS)
        u.ci.upp.ratioXY.NoS <- u.K.ratio.XY.NoS + (qnorm(0.975)*u.seK.ratio.XY.NoS)
        u.ci.low.ratioYX.NoS <- u.K.ratio.YX.NoS - (qnorm(0.975)*u.seK.ratio.YX.NoS)
        u.ci.upp.ratioYX.NoS <- u.K.ratio.YX.NoS + (qnorm(0.975)*u.seK.ratio.YX.NoS)
        
        u.ci.low.eggerXY.NoS <- u.K.egger.XY.NoS - (qnorm(0.975)*u.seK.egger.XY.NoS)
        u.ci.upp.eggerXY.NoS <- u.K.egger.XY.NoS + (qnorm(0.975)*u.seK.egger.XY.NoS)
        u.ci.low.eggerYX.NoS <- u.K.egger.YX.NoS - (qnorm(0.975)*u.seK.egger.YX.NoS)
        u.ci.upp.eggerYX.NoS <- u.K.egger.YX.NoS + (qnorm(0.975)*u.seK.egger.YX.NoS)
        
        #  DECISION RULES CD Ratio - .NoS unadjusted
        if((u.ci.low.ratioXY.NoS>0 & u.ci.upp.ratioXY.NoS<1)| (u.ci.low.ratioXY.NoS>(-1) & u.ci.upp.ratioXY.NoS<0)){
          if((u.ci.low.ratioYX.NoS>0 & u.ci.upp.ratioYX.NoS<1) | (u.ci.low.ratioYX.NoS>(-1) & u.ci.upp.ratioYX.NoS<0)){
            mat3[j,"biCDRatio.noS.u"] <- mat3[j,"biCDRatio.noS.u"]+1
            uR.noS <- 3
          }else{
            mat1[j,"biCDRatio.noS.u"] <- mat1[j,"biCDRatio.noS.u"]+1
            uR.noS <- 1
          }
        }else if((u.ci.low.ratioYX.NoS>0 & u.ci.upp.ratioYX.NoS<1) & (u.ci.low.ratioYX.NoS>(-1) & u.ci.upp.ratioYX.NoS<0)){
          if((u.ci.low.ratioXY.NoS>0 & u.ci.upp.ratioXY.NoS<1)| (u.ci.low.ratioXY.NoS>(-1) & u.ci.upp.ratioXY.NoS<0)){
            mat3[j,"biCDRatio.noS.u"] <- mat3[j,"biCDRatio.noS.u"]+1
            uR.noS <- 3
          }else{
            mat2[j,"biCDRatio.noS.u"] <- mat2[j,"biCDRatio.noS.u"]+1
            uR.noS <- 2
          }
        }else{
          mat3[j,"biCDRatio.noS.u"] <- mat3[j,"biCDRatio.noS.u"]+1
          uR.noS <- 3
        }
        
        #  DECISION RULES CD Egger - No S unadjusted
        if((u.ci.low.eggerXY.NoS>0 & u.ci.upp.eggerXY.NoS<1)| (u.ci.low.eggerXY.NoS>(-1) & u.ci.upp.eggerXY.NoS<0)){
          if((u.ci.low.eggerYX.NoS>0 & u.ci.upp.eggerYX.NoS<1) | (u.ci.low.eggerYX.NoS>(-1) & u.ci.upp.eggerYX.NoS<0)){
            mat3[j,"biCDEgger.noS.u"] <- mat3[j,"biCDEgger.noS.u"]+1
            uE.noS <- 3
          }else{
            mat1[j,"biCDEgger.noS.u"] <- mat1[j,"biCDEgger.noS.u"]+1
            uE.noS <- 1
          }
        }else if((u.ci.low.eggerYX.NoS>0 & u.ci.upp.eggerYX.NoS<1) & (u.ci.low.eggerYX.NoS>(-1) & u.ci.upp.eggerYX.NoS<0)){
          if((u.ci.low.eggerXY.NoS>0 & u.ci.upp.eggerXY.NoS<1)| (u.ci.low.eggerXY.NoS>(-1) & u.ci.upp.eggerXY.NoS<0)){
            mat3[j,"biCDEgger.noS.u"] <- mat3[j,"biCDEgger.noS.u"]+1
            uE.noS <- 3
          }else{
            mat2[j,"biCDEgger.noS.u"] <- mat2[j,"biCDEgger.noS.u"]+1
            uE.noS <- 2
          }
        }else{
          mat3[j,"biCDEgger.noS.u"] <- mat3[j,"biCDEgger.noS.u"]+1
          uE.noS <- 3
        }
        
        ########################################
        # Adjusted with screening
        ########################################
        aR.S <- NA
        aR.noS <- NA
        aE.S <- NA
        aE.noS <- NA
        
        a.K.ratio.XY.S <- res.biCD.A$CDRatio.XtoY.est.S
        a.seK.ratio.XY.S <- res.biCD.A$CDRatio.XtoY.se.S
        a.K.ratio.YX.S <- res.biCD.A$CDRatio.YtoX.est.S
        a.seK.ratio.YX.S <- res.biCD.A$CDRatio.YtoX.se.S
        
        a.K.egger.XY.S <- res.biCD.A$CDEgger.XtoY.est.S
        a.seK.egger.XY.S <- res.biCD.A$CDEgger.XtoY.se.S
        a.K.egger.YX.S <- res.biCD.A$CDEgger.YtoX.est.S
        a.seK.egger.YX.S <- res.biCD.A$CDEgger.YtoX.se.S
        
        a.ci.low.ratioXY.S <- a.K.ratio.XY.S - (qnorm(0.975)*a.seK.ratio.XY.S)
        a.ci.upp.ratioXY.S <- a.K.ratio.XY.S + (qnorm(0.975)*a.seK.ratio.XY.S)
        
        a.ci.low.ratioYX.S <- a.K.ratio.YX.S - (qnorm(0.975)*a.seK.ratio.YX.S)
        a.ci.upp.ratioYX.S <- a.K.ratio.YX.S + (qnorm(0.975)*a.seK.ratio.YX.S)
        
        a.ci.low.eggerXY.S <- a.K.egger.XY.S - (qnorm(0.975)*a.seK.egger.XY.S)
        a.ci.upp.eggerXY.S <- a.K.egger.XY.S + (qnorm(0.975)*a.seK.egger.XY.S)
        
        a.ci.low.eggerYX.S <- a.K.egger.YX.S - (qnorm(0.975)*a.seK.egger.YX.S)
        a.ci.upp.eggerYX.S <- a.K.egger.YX.S + (qnorm(0.975)*a.seK.egger.YX.S)
        
        #  DECISION RULES CD Ratio - .S - adjusted
        if((a.ci.low.ratioXY.S>0 & a.ci.upp.ratioXY.S<1)| (a.ci.low.ratioXY.S>-1 & a.ci.upp.ratioXY.S<0)){
          if((a.ci.low.ratioYX.S>0 & a.ci.upp.ratioYX.S<1) | (a.ci.low.ratioYX.S>-1 & a.ci.upp.ratioYX.S<0)){
            mat3[j,"biCDRatio.S.a"] <- mat3[j,"biCDRatio.S.a"]+1
            aR.S <- 3
          }else{
            mat1[j,"biCDRatio.S.a"] <- mat1[j,"biCDRatio.S.a"]+1
            aR.S <- 1
          }
        }else if((a.ci.low.ratioYX.S>0 & a.ci.upp.ratioYX.S<1) & (a.ci.low.ratioYX.S>-1 & a.ci.upp.ratioYX.S<0)){
          if((a.ci.low.ratioXY.S>0 & a.ci.upp.ratioXY.S<1)| (a.ci.low.ratioXY.S>-1 & a.ci.upp.ratioXY.S<0)){
            mat3[j,"biCDRatio.S.a"] <- mat3[j,"biCDRatio.S.a"]+1
            aR.S <- 3
          }else{
            mat2[j,"biCDRatio.S.a"] <- mat2[j,"biCDRatio.S.a"]+1
            aR.S <- 2
          }
        }else{
          mat3[j,"biCDRatio.S.a"] <- mat3[j,"biCDRatio.S.a"]+1
          aR.S <- 3
        }
        
        #  DECISION RULES CD Egger - adjusted
        if((a.ci.low.eggerXY.S>0 & a.ci.upp.eggerXY.S<1)| (a.ci.low.eggerXY.S>-1 & a.ci.upp.eggerXY.S<0)){
          if((a.ci.low.eggerYX.S>0 & a.ci.upp.eggerYX.S<1) | (a.ci.low.eggerYX.S>-1 & a.ci.upp.eggerYX.S<0)){
            mat3[j,"biCDEgger.S.a"] <- mat3[j,"biCDEgger.S.a"]+1
            aE.S <- 3
          }else{
            mat1[j,"biCDEgger.S.a"] <- mat1[j,"biCDEgger.S.a"]+1
            aE.S <- 1
          }
        }else if((a.ci.low.eggerYX.S>0 & a.ci.upp.eggerYX.S<1) & (a.ci.low.eggerYX.S>-1 & a.ci.upp.eggerYX.S<0)){
          if((a.ci.low.eggerXY.S>0 & a.ci.upp.eggerXY.S<1)| (a.ci.low.eggerXY.S>-1 & a.ci.upp.eggerXY.S<0)){
            mat3[j,"biCDEgger.S.a"] <- mat3[j,"biCDEgger.S.a"]+1
            aE.S <- 3
          }else{
            mat2[j,"biCDEgger.S.a"] <- mat2[j,"biCDEgger.S.a"]+1
            aE.S <- 2
          }
        }else{
          mat3[j,"biCDEgger.S.a"] <- mat3[j,"biCDEgger.S.a"]+1
          aE.S <- 3
        }
        
        # save K & SE for no screening
        a.K.ratio.XY.NoS <- res.biCD.A$CDRatio.XtoY.est.NoS
        a.seK.ratio.XY.NoS <- res.biCD.A$CDRatio.XtoY.se.NoS
        a.K.ratio.YX.NoS <- res.biCD.A$CDRatio.YtoX.est.NoS
        a.seK.ratio.YX.NoS <- res.biCD.A$CDRatio.YtoX.se.NoS
        
        a.K.egger.XY.NoS <- res.biCD.A$CDEgger.XtoY.est.NoS
        a.seK.egger.XY.NoS <- res.biCD.A$CDEgger.XtoY.se.NoS
        a.K.egger.YX.NoS <- res.biCD.A$CDEgger.YtoX.est.NoS
        a.seK.egger.YX.NoS <- res.biCD.A$CDEgger.YtoX.se.NoS
        
        a.ci.low.ratioXY.NoS <- a.K.ratio.XY.NoS - (qnorm(0.975)*a.seK.ratio.XY.NoS)
        a.ci.upp.ratioXY.NoS <- a.K.ratio.XY.NoS + (qnorm(0.975)*a.seK.ratio.XY.NoS)
        
        a.ci.low.ratioYX.NoS <- a.K.ratio.YX.NoS - (qnorm(0.975)*a.seK.ratio.YX.NoS)
        a.ci.upp.ratioYX.NoS <- a.K.ratio.YX.NoS + (qnorm(0.975)*a.seK.ratio.YX.NoS)
        
        a.ci.low.eggerXY.NoS <- a.K.egger.XY.NoS - (qnorm(0.975)*a.seK.egger.XY.NoS)
        a.ci.upp.eggerXY.NoS <- a.K.egger.XY.NoS + (qnorm(0.975)*a.seK.egger.XY.NoS)
        
        a.ci.low.eggerYX.NoS <- a.K.egger.YX.NoS - (qnorm(0.975)*a.seK.egger.YX.NoS)
        a.ci.upp.eggerYX.NoS <- a.K.egger.YX.NoS + (qnorm(0.975)*a.seK.egger.YX.NoS)
        
        #  DECISION RULES CD Ratio - .NoS
        if((a.ci.low.ratioXY.NoS>0 & a.ci.upp.ratioXY.NoS<1)| (a.ci.low.ratioXY.NoS>-1 & a.ci.upp.ratioXY.NoS<0)){
          if((a.ci.low.ratioYX.NoS>0 & a.ci.upp.ratioYX.NoS<1) | (a.ci.low.ratioYX.NoS>-1 & a.ci.upp.ratioYX.NoS<0)){
            mat3[j,"biCDRatio.noS.a"] <- mat3[j,"biCDRatio.noS.a"]+1
            aR.noS <- 3
          }else{
            mat1[j,"biCDRatio.noS.a"] <- mat1[j,"biCDRatio.noS.a"]+1
            aR.noS <- 1
          }
        }else if((a.ci.low.ratioYX.NoS>0 & a.ci.upp.ratioYX.NoS<1) & (a.ci.low.ratioYX.NoS>-1 & a.ci.upp.ratioYX.NoS<0)){
          if((a.ci.low.ratioXY.NoS>0 & a.ci.upp.ratioXY.NoS<1)| (a.ci.low.ratioXY.NoS>-1 & a.ci.upp.ratioXY.NoS<0)){
            mat3[j,"biCDRatio.noS.a"] <- mat3[j,"biCDRatio.noS.a"]+1
            aR.noS <- 3
          }else{
            mat2[j,"biCDRatio.noS.a"] <- mat2[j,"biCDRatio.noS.a"]+1
            aR.noS <- 2
          }
        }else{
          mat3[j,"biCDRatio.noS.a"] <- mat3[j,"biCDRatio.noS.a"]+1
          aR.noS <- 3
        }
        
        #  DECISION RULES CD Egger
        if((a.ci.low.eggerXY.NoS>0 & a.ci.upp.eggerXY.NoS<1)| (a.ci.low.eggerXY.NoS>-1 & a.ci.upp.eggerXY.NoS<0)){
          if((a.ci.low.eggerYX.NoS>0 & a.ci.upp.eggerYX.NoS<1) | (a.ci.low.eggerYX.NoS>-1 & a.ci.upp.eggerYX.NoS<0)){
            mat3[j,"biCDEgger.noS.a"] <- mat3[j,"biCDEgger.noS.a"]+1
            aE.noS <- 3
          }else{
            mat1[j,"biCDEgger.noS.a"] <- mat1[j,"biCDEgger.noS.a"]+1
            aE.noS <- 1
          }
        }else if((a.ci.low.eggerYX.NoS>0 & a.ci.upp.eggerYX.NoS<1) & (a.ci.low.eggerYX.NoS>-1 & a.ci.upp.eggerYX.NoS<0)){
          if((a.ci.low.eggerXY.NoS>0 & a.ci.upp.eggerXY.NoS<1)| (a.ci.low.eggerXY.NoS>-1 & a.ci.upp.eggerXY.NoS<0)){
            mat3[j,"biCDEgger.noS.a"] <- mat3[j,"biCDEgger.noS.a"]+1
            aE.noS <- 3
          }else{
            mat2[j,"biCDEgger.noS.a"] <- mat2[j,"biCDEgger.noS.a"]+1
            aE.noS <- 2
          }
        }else{
          mat3[j,"biCDEgger.noS.a"] <- mat3[j,"biCDEgger.noS.a"]+1
          aE.noS <- 3
        }
        
        ########################################################
        # fill in matDiff for CD s and no s approaches
        ########################################################
        matDiff[j,"biCDRatioXY.S"] <- matDiff[j,"biCDRatioXY.S"] + abs(a.K.ratio.XY.S - u.K.ratio.XY.S)
        matDiff[j,"biCDRatioYX.S"] <- matDiff[j,"biCDRatioYX.S"] + abs(a.K.ratio.YX.S - u.K.ratio.YX.S)
        
        matDiff[j,"biCDRatioXY.noS"] <- matDiff[j,"biCDRatioXY.noS"] + abs(a.K.ratio.XY.NoS - u.K.ratio.XY.NoS)
        matDiff[j,"biCDRatioYX.noS"] <- matDiff[j,"biCDRatioYX.noS"] + abs(a.K.ratio.YX.NoS - u.K.ratio.YX.NoS)
        
        matDiff[j,"biCDEggerXY.S"] <- matDiff[j,"biCDEggerXY.S"] + abs(a.K.egger.XY.S - u.K.egger.XY.S)
        matDiff[j,"biCDEggerYX.S"] <- matDiff[j,"biCDEggerYX.S"] + abs(a.K.egger.YX.S - u.K.egger.YX.S)
        
        matDiff[j,"biCDEggerXY.noS"] <- matDiff[j,"biCDEggerXY.noS"] + abs(a.K.egger.XY.NoS - u.K.egger.XY.NoS)
        matDiff[j,"biCDEggerYX.noS"] <- matDiff[j,"biCDEggerYX.noS"] + abs(a.K.egger.YX.NoS - u.K.egger.YX.NoS)
        
        ##################################################################
        # Fill in matSame for Bi CD Ratio, Bi CD Egger
        ##################################################################
        if(uR.S==aR.S){matSame[j,"biCDRatio.S"] <- matSame[j,"biCDRatio.S"] + 1}
        
        if(uR.noS==aR.noS){matSame[j,"biCDRatio.noS"] <- matSame[j,"biCDRatio.noS"] + 1}
        
        if(uE.S==aE.S){matSame[j,"biCDEgger.S"] <- matSame[j,"biCDEgger.S"] + 1}
        
        if(uE.noS==aE.noS){matSame[j,"biCDEgger.noS"] <- matSame[j,"biCDEgger.noS"]+1}
        
        ##################################
        # CD-Ratio, CD-Egger, CD-GLS
        ##################################
        # betas for CD approaches
        # GX --> X --> Y
        u.betaGX_X <- NULL
        u.seGX_X <- NULL
        u.betaGX_Y <- NULL
        u.seGX_Y <- NULL
        for(ll in 1:ncol(matGX)){
          u.betaGX_X <- c(u.betaGX_X, summary(lm(x~matGX[,ll]))$coef[2,1])
          u.seGX_X <- c(u.seGX_X, summary(lm(x~matGX[,ll]))$coef[2,2])
          
          u.betaGX_Y <- c(u.betaGX_Y, summary(lm(y~matGX[,ll]))$coef[2,1])
          u.seGX_Y <- c(u.seGX_Y, summary(lm(y~matGX[,ll]))$coef[2,2])
        }
        
        a.betaGX_X <- NULL
        a.seGX_X <- NULL
        a.betaGX_Y <- NULL
        a.seGX_Y <- NULL
        for(ll in 1:ncol(matGX)){
          a.betaGX_X <- c(a.betaGX_X, summary(lm(x~matGX[,ll]+matCX))$coef[2,1])
          a.seGX_X <- c(a.seGX_X, summary(lm(x~matGX[,ll]+matCX))$coef[2,2])
          
          a.betaGX_Y <- c(a.betaGX_Y, summary(lm(y~matGX[,ll]+matCY))$coef[2,1])
          a.seGX_Y <- c(a.seGX_Y, summary(lm(y~matGX[,ll]+matCY))$coef[2,2])
        }
        
        ######################################################################################################
        # Unadjusted -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        uR <- NA
        uG <- NA
        uE <- NA
        
        u.tmp.mat <- matrix(NA,ncol=12,nrow=length(u.betaGX_X))
        u.tmp.df <- as.data.frame(u.tmp.mat)
        colnames(u.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        u.tmp.df$rsid<-paste0("SNP", c(1:length(u.betaGX_X)))
        u.tmp.df$beta_T1 <- u.betaGX_X
        u.tmp.df$se_T1 <- u.seGX_X
        u.tmp.df$N_T1 <- n
        u.tmp.df$beta_T2 <- u.betaGX_Y
        u.tmp.df$se_T2 <- u.seGX_Y
        u.tmp.df$N_T2 <- n
        
        u.pruned1 <- list(sig_part=u.tmp.df)
        u.cd3 <- CD_3_methods_Independent(u.pruned1$sig_part)
        
        u.ratio.YX <- u.cd3$CD_Ratio_result$T1toT2
        u.ratio.XY <- u.cd3$CD_Ratio_result$T2toT1
        u.egger.YX<- u.cd3$CD_Egger_result$T1toT2
        u.egger.XY <- u.cd3$CD_Egger_result$T2toT1
        u.gls.YX <- u.cd3$CD_GLS_result$T1toT2
        u.gls.XY <- u.cd3$CD_GLS_result$T2toT1
        
        # Confidence intervals for K - needed for use in decision rules
        # CD-ratio CIs
        u.lowerCIyx.cdRatio <- u.ratio.YX["K"] - qnorm(1-sig.level/2)*u.ratio.YX["se(K)"]
        u.upperCIyx.cdRatio <- u.ratio.YX["K"] + qnorm(1-sig.level/2)*u.ratio.YX["se(K)"]
        u.lowerCIxy.cdRatio <- u.ratio.XY["K"] - qnorm(1-sig.level/2)*u.ratio.XY["se(K)"]
        u.upperCIxy.cdRatio <- u.ratio.XY["K"] + qnorm(1-sig.level/2)*u.ratio.XY["se(K)"]
        
        # CD-Egger CIs
        u.lowerCIyx.cdEgger <- u.egger.YX["K"] - qnorm(1-sig.level/2)*u.egger.YX["se(K)"]
        u.upperCIyx.cdEgger <- u.egger.YX["K"] + qnorm(1-sig.level/2)*u.egger.YX["se(K)"]
        u.lowerCIxy.cdEgger <- u.egger.XY["K"] - qnorm(1-sig.level/2)*u.egger.XY["se(K)"]
        u.upperCIxy.cdEgger <- u.egger.XY["K"] + qnorm(1-sig.level/2)*u.egger.XY["se(K)"]
        
        # CD-GLS CIs
        u.lowerCIyx.cdGLS <- u.gls.YX["K"] - qnorm(1-sig.level/2)*u.gls.YX["se(K)"]
        u.upperCIyx.cdGLS <- u.gls.YX["K"] + qnorm(1-sig.level/2)*u.gls.YX["se(K)"]
        u.lowerCIxy.cdGLS <- u.gls.XY["K"] - qnorm(1-sig.level/2)*u.gls.XY["se(K)"]
        u.upperCIxy.cdGLS <- u.gls.XY["K"] + qnorm(1-sig.level/2)*u.gls.XY["se(K)"]
        
        # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
        # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
        # case 3: neither otherwise
        # CD-Ratio decisions
        if(u.upperCIxy.cdRatio < u.lowerCIxy.cdRatio){
          stop("LowerCIxy > UpperCIxy")}
        if(u.upperCIyx.cdRatio < u.lowerCIyx.cdRatio){
          stop("LowerCIyx > UpperCIyx")}
        
        if(u.upperCIxy.cdEgger < u.lowerCIxy.cdEgger){
          stop("LowerCIxyEgger > UpperCIxyEgger")}
        if(u.upperCIyx.cdEgger < u.lowerCIyx.cdEgger){
          stop("LowerCIyxEgger > UpperCIyxEgger")}
        
        if(u.upperCIxy.cdGLS < u.lowerCIxy.cdGLS){
          stop("LowerCIxy > UpperCIxy")}
        if(u.upperCIyx.cdGLS < u.lowerCIyx.cdGLS){
          stop("LowerCIyx > UpperCIyx")}
        
        # CD-Ratio decisions
        if((u.upperCIxy.cdRatio < (-1) | u.lowerCIxy.cdRatio>1) &
           ((u.lowerCIyx.cdRatio>=(-1) & u.upperCIyx.cdRatio<0) | (u.lowerCIyx.cdRatio>(0) & u.upperCIyx.cdRatio<=1)) ){
          mat1[j, "CDRatio.u"] <- mat1[j, "CDRatio.u"] + 1
          uR <- 1
        }else if((u.upperCIyx.cdRatio < (-1) | u.lowerCIyx.cdRatio>1) &
                 ((u.lowerCIxy.cdRatio>=(-1) & u.upperCIxy.cdRatio< 0) | (u.lowerCIxy.cdRatio>(0) & u.upperCIxy.cdRatio<= 1)) ){
          mat2[j, "CDRatio.u"] <- mat2[j, "CDRatio.u"] + 1
          uR <- 2
        }else{
          mat3[j, "CDRatio.u"] <- mat3[j, "CDRatio.u"] + 1
          uR <- 3
        }
        
        # CD-Egger decisions
        if((u.upperCIxy.cdEgger < (-1)  | u.lowerCIxy.cdEgger>1) &
           ((u.lowerCIyx.cdEgger>=(-1) & u.upperCIyx.cdEgger<0) | (u.lowerCIyx.cdEgger>(0) & u.upperCIyx.cdEgger<=1))){
          mat1[j, "CDEgger.u"] <- mat1[j, "CDEgger.u"] + 1
          uE <- 1
        }else if((u.upperCIyx.cdEgger < (-1) | u.lowerCIyx.cdEgger>1) &
                 ((u.lowerCIxy.cdEgger>=(-1)  & u.upperCIxy.cdEgger< 0) | (u.lowerCIxy.cdEgger>(0)  & u.upperCIxy.cdEgger<=1))){
          mat2[j, "CDEgger.u"] <- mat2[j, "CDEgger.u"] + 1
          uE <- 2
        }else{
          mat3[j, "CDEgger.u"] <- mat3[j, "CDEgger.u"] + 1
          uE <- 3
        }
        
        # CD-GLS decisions
        if((u.upperCIxy.cdGLS < (-1) | u.lowerCIxy.cdGLS>1) &
           ((u.lowerCIyx.cdGLS>= (-1) & u.upperCIyx.cdGLS<0) | (u.lowerCIyx.cdGLS> (0) & u.upperCIyx.cdGLS<=1)) ){
          mat1[j, "CDGLS.u"] <- mat1[j, "CDGLS.u"] + 1
          uG <- 1
        }else if((u.upperCIyx.cdGLS < (-1) | u.lowerCIyx.cdGLS>1) &
                 ((u.lowerCIxy.cdGLS>= (-1) & u.upperCIxy.cdGLS<0) | (u.lowerCIxy.cdGLS> (0) & u.upperCIxy.cdGLS<=1)) ){
          mat2[j, "CDGLS.u"] <- mat2[j, "CDGLS.u"] + 1
          uG <- 2
        }else{
          mat3[j, "CDGLS.u"] <- mat3[j, "CDGLS.u"] + 1
          uG <- 3
        }
        
        
        ######################################################################################################
        # Adjusted -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        aR <- NA
        aG <- NA
        aE <- NA
        
        a.tmp.mat <- matrix(NA,ncol=12,nrow=length(a.betaGX_X))
        a.tmp.df <- as.data.frame(a.tmp.mat)
        colnames(a.tmp.df) <- c("chr", "pos", "rsid", "A1", "A2", "beta_T1",
                                "se_T1", "N_T1", "beta_T2", "se_T2", "N_T2", "loci")
        a.tmp.df$rsid<-paste0("SNP", c(1:length(a.betaGX_X)))
        a.tmp.df$beta_T1 <- a.betaGX_X
        a.tmp.df$se_T1 <- a.seGX_X
        a.tmp.df$N_T1 <- n
        a.tmp.df$beta_T2 <- a.betaGX_Y
        a.tmp.df$se_T2 <- a.seGX_Y
        a.tmp.df$N_T2 <- n
        
        a.pruned1 <- list(sig_part=a.tmp.df)
        a.cd3 <- CD_3_methods_Independent(a.pruned1$sig_part)
        
        a.ratio.YX <- a.cd3$CD_Ratio_result$T1toT2
        a.ratio.XY <- a.cd3$CD_Ratio_result$T2toT1
        a.egger.YX<- a.cd3$CD_Egger_result$T1toT2
        a.egger.XY <- a.cd3$CD_Egger_result$T2toT1
        a.gls.YX <- a.cd3$CD_GLS_result$T1toT2
        a.gls.XY <- a.cd3$CD_GLS_result$T2toT1
        
        # Confidence intervals for K - needed for use in decision rules
        # CD-ratio CIs
        a.lowerCIyx.cdRatio <- a.ratio.YX["K"] - qnorm(1-sig.level/2)*a.ratio.YX["se(K)"]
        a.upperCIyx.cdRatio <- a.ratio.YX["K"] + qnorm(1-sig.level/2)*a.ratio.YX["se(K)"]
        a.lowerCIxy.cdRatio <- a.ratio.XY["K"] - qnorm(1-sig.level/2)*a.ratio.XY["se(K)"]
        a.upperCIxy.cdRatio <- a.ratio.XY["K"] + qnorm(1-sig.level/2)*a.ratio.XY["se(K)"]
        
        # CD-Egger CIs
        a.lowerCIyx.cdEgger <- a.egger.YX["K"] - qnorm(1-sig.level/2)*a.egger.YX["se(K)"]
        a.upperCIyx.cdEgger <- a.egger.YX["K"] + qnorm(1-sig.level/2)*a.egger.YX["se(K)"]
        a.lowerCIxy.cdEgger <- a.egger.XY["K"] - qnorm(1-sig.level/2)*a.egger.XY["se(K)"]
        a.upperCIxy.cdEgger <- a.egger.XY["K"] + qnorm(1-sig.level/2)*a.egger.XY["se(K)"]
        
        # CD-GLS CIs
        a.lowerCIyx.cdGLS <- a.gls.YX["K"] - qnorm(1-sig.level/2)*a.gls.YX["se(K)"]
        a.upperCIyx.cdGLS <- a.gls.YX["K"] + qnorm(1-sig.level/2)*a.gls.YX["se(K)"]
        a.lowerCIxy.cdGLS <- a.gls.XY["K"] - qnorm(1-sig.level/2)*a.gls.XY["se(K)"]
        a.upperCIxy.cdGLS <- a.gls.XY["K"] + qnorm(1-sig.level/2)*a.gls.XY["se(K)"]
        
        # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
        # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
        # case 3: neither otherwise
        # CD-Ratio decisions
        if(a.upperCIxy.cdRatio < a.lowerCIxy.cdRatio){
          stop("LowerCIxy > UpperCIxy")}
        if(a.upperCIyx.cdRatio < a.lowerCIyx.cdRatio){
          stop("LowerCIyx > UpperCIyx")}
        
        if(a.upperCIxy.cdEgger < a.lowerCIxy.cdEgger){
          stop("LowerCIxyEgger > UpperCIxyEgger")}
        if(a.upperCIyx.cdEgger < a.lowerCIyx.cdEgger){
          stop("LowerCIyxEgger > UpperCIyxEgger")}
        
        if(a.upperCIxy.cdGLS < a.lowerCIxy.cdGLS){
          stop("LowerCIxy > UpperCIxy")}
        if(a.upperCIyx.cdGLS < a.lowerCIyx.cdGLS){
          stop("LowerCIyx > UpperCIyx")}
        
        # CD-Ratio decisions
        if((a.upperCIxy.cdRatio < (-1) | a.lowerCIxy.cdRatio>1) &
           ((a.lowerCIyx.cdRatio>=(-1) & a.upperCIyx.cdRatio<0) | (a.lowerCIyx.cdRatio>(0) & a.upperCIyx.cdRatio<=1)) ){
          mat1[j, "CDRatio.a"] <- mat1[j, "CDRatio.a"] + 1
          aR <- 1
        }else if((a.upperCIyx.cdRatio < (-1) | a.lowerCIyx.cdRatio>1) &
                 ((a.lowerCIxy.cdRatio>=(-1) & a.upperCIxy.cdRatio< 0) | (a.lowerCIxy.cdRatio>(0) & a.upperCIxy.cdRatio<= 1)) ){
          mat2[j, "CDRatio.a"] <- mat2[j, "CDRatio.a"] + 1
          aR <- 2
        }else{
          mat3[j, "CDRatio.a"] <- mat3[j, "CDRatio.a"] + 1
          aR <- 3
        }
        
        # CD-Egger decisions
        if((a.upperCIxy.cdEgger < (-1)  | a.lowerCIxy.cdEgger>1) &
           ((a.lowerCIyx.cdEgger>=(-1) & a.upperCIyx.cdEgger<0) | (a.lowerCIyx.cdEgger>(0) & a.upperCIyx.cdEgger<=1))){
          mat1[j, "CDEgger.a"] <- mat1[j, "CDEgger.a"] + 1
          aE <- 1
        }else if((a.upperCIyx.cdEgger < (-1) | a.lowerCIyx.cdEgger>1) &
                 ((a.lowerCIxy.cdEgger>=(-1)  & a.upperCIxy.cdEgger< 0) | (a.lowerCIxy.cdEgger>(0)  & a.upperCIxy.cdEgger<=1))){
          mat2[j, "CDEgger.a"] <- mat2[j, "CDEgger.a"] + 1
          aE <- 2
        }else{
          mat3[j, "CDEgger.a"] <- mat3[j, "CDEgger.a"] + 1
          aE <- 3
        }
        
        # CD-GLS decisions
        if((a.upperCIxy.cdGLS < (-1) | a.lowerCIxy.cdGLS>1) &
           ((a.lowerCIyx.cdGLS>= (-1) & a.upperCIyx.cdGLS<0) | (a.lowerCIyx.cdGLS> (0) & a.upperCIyx.cdGLS<=1)) ){
          mat1[j, "CDGLS.a"] <- mat1[j, "CDGLS.a"] + 1
          aG <- 1
        }else if((a.upperCIyx.cdGLS < (-1) | a.lowerCIyx.cdGLS>1) &
                 ((a.lowerCIxy.cdGLS>= (-1) & a.upperCIxy.cdGLS<0) | (a.lowerCIxy.cdGLS> (0) & a.upperCIxy.cdGLS<=1)) ){
          mat2[j, "CDGLS.a"] <- mat2[j, "CDGLS.a"] + 1
          aG <- 2
        }else{
          mat3[j, "CDGLS.a"] <- mat3[j, "CDGLS.a"] + 1
          aG <- 3
        }
        
        
        ########################################################
        # fill in matDiff for CD s and no s approaches
        ########################################################
        matDiff[j,"CDRatioXY"] <- matDiff[j,"CDRatioXY"] + abs(a.ratio.XY["K"] - u.ratio.XY["K"] )
        matDiff[j,"CDRatioYX"] <- matDiff[j,"CDRatioYX"] + abs(a.ratio.YX["K"] - u.ratio.YX["K"] )
        
        matDiff[j,"CDEggerXY"] <- matDiff[j,"CDEggerXY"] + abs(a.egger.XY["K"] - u.egger.XY["K"])
        matDiff[j,"CDEggerYX"] <- matDiff[j,"CDEggerYX"] + abs(a.egger.YX["K"] - u.egger.YX["K"])
        
        matDiff[j,"CDGLSXY"] <- matDiff[j,"CDGLSXY"] + abs(a.gls.XY["K"] - u.gls.XY["K"])
        matDiff[j,"CDGLSYX"] <- matDiff[j,"CDGLSYX"] + abs(a.gls.YX["K"] - u.gls.YX["K"])
        
        
        ##################################################################
        # Fill in matSame for CD s and no s approaches
        ##################################################################
        if(uR==aR){matSame[j,"CDRatio"] <- matSame[j,"CDRatio"] + 1}
        
        if(uG==aG){matSame[j,"CDGLS"] <- matSame[j,"CDGLS"] + 1}
        
        if(uE==aE){matSame[j,"CDEgger"] <- matSame[j,"CDEgger"] + 1}
        
      }# end cov loop
    }# end sims loop
    
    
    mat_total1 <- cbind(mat1[,1], mat1[,-1]/nSims)
    mat_total2 <- cbind(mat2[,1], mat2[,-1]/nSims)
    mat_total3 <- cbind(mat3[,1], mat3[,-1]/nSims)
    mat_Diff <- cbind(matDiff[,1], matDiff[,-1]/nSims)
    mat_Same <- cbind(matSame[,1], matSame[,-1]/nSims)
    
    
    
    if(plot.pdf){
      
      #########################################
      # Difference in K plot - X --> Y
      #########################################
      pdf(paste(plot.name,"_Diff_XY.pdf", sep = ""))
      plot(-2,-2,xlim=c(min(deltaC),max(deltaC)),ylim=c(0,max(mat_Diff[,2:ncol(mat_Diff)])),main="",xlab=expression(delta[C]),ylab="Absolute difference")
      
      ###########
      lines(deltaC[,1],mat_Diff[,"CDcML.S.DP.XY"],col="steelblue1",pch=3,lty=6,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDcML.S.noDP.XY"],col="darkorchid4",pch=1,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDcML.noS.DP.XY"],col="dodgerblue",pch=4,lty=3,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDcML.noS.noDP.XY"],col="mediumorchid1",pch=6,lty=1,type="b", lwd=1.2)
      ###########      
      lines(deltaC[,1],mat_Diff[,"MRcML.S.DP.XY"],col="deeppink2",pch=1,lty=1,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"MRcML.S.noDP.XY"],col="firebrick1",pch=3,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"MRcML.noS.DP.XY"],col="darkorange1",pch=4,lty=3,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"MRcML.noS.noDP.XY"],col="orange1",pch=0,lty=6,type="b", lwd=1.2)
      #########
      lines(deltaC[,1],mat_Diff[,"biCDRatioXY.S"],col="goldenrod3",pch=2,lty=1,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"biCDRatioXY.noS"],col="goldenrod1",pch=3,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"biCDEggerXY.S"],col="springgreen4",pch=0,lty=3,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"biCDEggerXY.noS"],col="seagreen4",pch=6,lty=6,type="b", lwd=1.2)
      #########
      lines(deltaC[,1],mat_Diff[,"CDRatioXY"],col="slategray1",pch=2,lty=1,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDEggerXY"],col="grey52",pch=3,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDGLSXY"],col="black",pch=6,lty=6,type="b", lwd=1.2)
      
      dev.off()      
      
      
      
      #########################################
      # Difference in K plot - Y --> X
      #########################################
      pdf(paste(plot.name,"_Diff_YX.pdf", sep = ""))
      plot(-2,-2,xlim=c(min(deltaC),max(deltaC)),ylim=c(0,max(mat_Diff[,2:ncol(mat_Diff)])),main="",xlab=expression(delta[C]),ylab="Absolute difference")
      
      ###########
      lines(deltaC[,1],mat_Diff[,"CDcML.S.DP.YX"],col="steelblue1",pch=3,lty=6,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDcML.S.noDP.YX"],col="darkorchid4",pch=1,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDcML.noS.DP.YX"],col="dodgerblue",pch=4,lty=3,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDcML.noS.noDP.YX"],col="mediumorchid1",pch=6,lty=1,type="b", lwd=1.2)
      ##########
      lines(deltaC[,1],mat_Diff[,"MRcML.S.DP.YX"],col="deeppink2",pch=1,lty=1,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"MRcML.S.noDP.YX"],col="firebrick1",pch=3,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"MRcML.noS.DP.YX"],col="darkorange1",pch=4,lty=3,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"MRcML.noS.noDP.YX"],col="orange1",pch=0,lty=6,type="b", lwd=1.2)
      #########
      lines(deltaC[,1],mat_Diff[,"biCDRatioYX.S"],col="goldenrod3",pch=2,lty=1,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"biCDRatioYX.noS"],col="goldenrod1",pch=3,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"biCDEggerYX.S"],col="springgreen4",pch=0,lty=3,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"biCDEggerYX.noS"],col="seagreen4",pch=6,lty=6,type="b", lwd=1.2)
      #########
      lines(deltaC[,1],mat_Diff[,"CDRatioYX"],col="slategray1",pch=2,lty=1,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDEggerYX"],col="grey52",pch=3,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Diff[,"CDGLSYX"],col="black",pch=6,lty=6,type="b", lwd=1.2)
      
      dev.off()          
      
      #########################################
      # Proportion the same plot
      #########################################
      pdf(paste(plot.name,"_Same.pdf", sep = ""))
      plot(-2,-2,xlim=c(min(deltaC),max(deltaC)),ylim=c(0,1.1),main="",xlab=expression(beta[C]),ylab="Proportion of Simulations")
      
      ###########
      lines(deltaC[,1],mat_Same[,"CDcML.S.DP"],col="steelblue1",pch=3,lty=6,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"CDcML.S.noDP"],col="darkorchid4",pch=1,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"CDcML.noS.DP"],col="dodgerblue",pch=4,lty=3,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"CDcML.noS.noDP"],col="mediumorchid1",pch=2,lty=1,type="b", lwd=1.2)
      
      lines(deltaC[,1],mat_Same[,"MRcML.S.DP"],col="deeppink2",pch=1,lty=1,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"MRcML.S.noDP"],col="firebrick1",pch=3,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"MRcML.noS.DP"],col="darkorange1",pch=4,lty=3,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"MRcML.noS.noDP"],col="orange1",pch=0,lty=6,type="b", lwd=1.2)
      
      #########
      lines(deltaC[,1],mat_Same[,"biCDRatio.S"],col="goldenrod3",pch=2,lty=1,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"biCDRatio.noS"],col="goldenrod1",pch=3,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"biCDEgger.S"],col="springgreen4",pch=0,lty=3,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"biCDEgger.noS"],col="seagreen4",pch=6,lty=6,type="b", lwd=1.2)
      #########
      lines(deltaC[,1],mat_Same[,"CDRatio"],col="slategray1",pch=2,lty=1,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"CDEgger"],col="grey52",pch=3,lty=2,type="b", lwd=1.2)
      lines(deltaC[,1],mat_Same[,"CDGLS"],col="black",pch=6,lty=6,type="b", lwd=1.2)
      
      dev.off()      
      
      
      
    } # end plot if statement
    if(table.sv){
      write.table(matSame, file=paste0(plot.name,"_matSame.txt"), quote=F, row.names=F)
      write.table(matDiff, file=paste0(plot.name,"_matDiff.txt"), quote=F, row.names=F)
      write.table(mat1, file=paste0(plot.name,"_mat1.txt"), quote=F, row.names=F)
      write.table(mat2, file=paste0(plot.name,"_mat2.txt"), quote=F, row.names=F)
      write.table(mat3, file=paste0(plot.name,"_mat3.txt"), quote=F, row.names=F)
    }
    
    
    return(list("mat1"=mat1, "mat2"=mat2, "mat3"=mat3, "matDiff"=matDiff, "matSame"= matSame))
    
    
  }
