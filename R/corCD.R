corCD <-
function(n=100000, nSNPx = 10, MAFx = c(rep(0.5,  10)), nSNPy = 10, MAFy = c(rep(0.5, 10)),
           nCov = 3, meanC = rep(0, 3), sdC = c(rep(1, 3)), deltaC =  c(rep(0.005,3)), betaC=c(rep(0.005,3)),
           betaGY = c(rep(0.2, 10)), betaGX = c(rep(0, 10)), deltaGX = c(rep(0.2, 10)), betaX=0.2,
           sdX=1,sdY=1, SEED=1, sig.level=0.05, nSims=5,  plot.pdf=FALSE, plot.name="plot", table.sv=FALSE, 
           approaches=c("All")){
    
    if(length(MAFx)!=nSNPx){stop("length(MAFx) must equal nSNPx.")}
    if(length(MAFy)!=nSNPy){stop("length(MAFy) must equal nSNPy")}
    if(length(deltaGX)!=nSNPx){stop("length(deltaGX) must equal nSNPx")}
    if(length(betaGY)!=nSNPy){stop("length(betaGY) must equal nSNPy")}
    if(length(betaGX)!=nSNPy){stop("length(betaGX) must equal nSNPy")}
    
    if(length(meanC)!=nCov){stop("length(meanC) must equal nCov")}
    if(length(sdC)!=nCov){stop("length(sdC) must equal nCov")}
    if(length(deltaC)!=nCov){stop("ncol(deltaC) must equal nCov")}
    if(length(betaC)!=nCov){stop("length(betaC) must equal nCov")}
    
    # load libraries
    library(MRCD)
    library(MRcML)
    library(BiDirectCausal)
    library(RColorBrewer)

    ################################################################################
    # Matrix to save Results
    ################################################################################
    #save results for type 1 error rate betaX=0 and power betaX>0
    
    All <- c("CDcML", "MRcML", "CDRatio", "CDEgger", "CDGLS")
    
    colnames.1 <- NULL
    if("All" %in% approaches){colnames.1 <- c(colnames.1, All)}
    if("CDcML" %in% approaches){colnames.1 <- c(colnames.1, "CDcML")}
    if("MRcML" %in% approaches){colnames.1 <- c(colnames.1, "MRcML")}
    
    if("CDRatio" %in% approaches){colnames.1 <- c(colnames.1, "CDRatio")}
    if("CDEgger" %in% approaches){colnames.1 <- c(colnames.1, "CDEgger")}
    if("CDGLS" %in% approaches){colnames.1 <- c(colnames.1, "CDGLS")}
    
    colnames.mat1 <- unique(colnames.1)
    
    
    matR <- matrix(0, ncol=length(colnames.mat1), nrow=9)
    
    colnames(matR) <- colnames.mat1
    rownames(matR) <- c("slr1_mlr1", "slr1_mlr2", "slr1_mlr3", "slr2_mlr1", "slr2_mlr2", "slr2_mlr3", "slr3_mlr1", "slr3_mlr2", "slr3_mlr3")
    
    #storing matrices/vectors
    corGX<-rep(0,nSims)
    corGY<-rep(0,nSims)
    
    for(i in 1:nSims){
      printCut<-1
      if(floor(i/printCut)==ceiling(i/printCut)){print(paste(i,"of",nSims, "simulations"))}
      set.seed(SEED+i-1)
      
      ###################################
      # matrix for deltaGX, deltaC, betaGX, betaGY, betaC
      deltaGX<-matrix(deltaGX, ncol=1, nrow=nSNPx)
      betaC<-matrix(betaC, ncol=1, nrow=nCov)
      betaGX<- matrix(betaGX,ncol=1,nrow=nSNPy)
      betaGY<- matrix(betaGY,ncol=1,nrow=nSNPy)
      deltaC<-matrix(deltaC, ncol=1, nrow=nCov)
      
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
      
      # ONLY ONE SET OF COV - remove
      # create Covariate matrix
      matC<-matrix(0, ncol=nCov,nrow=n)
      for(cc in 1:nCov){
        matC[,cc]<-rnorm(n, meanC[cc], sdC[cc])
      }
      
      
      x<-rnorm(n, (matGX%*%deltaGX + matC%*%deltaC), sdX) # trait X
      y<- rnorm(n, (matGY%*%betaGY + matGX%*%betaGX+ matC%*%betaC+betaX*x), sdY) # trait X 
      
      
      
      ###################################
      # correlations
      ###################################

      # simple linear regression (SLR)
      modelXS<-summary(lm(x~matGX[,1]))$coef[2,c(1,2)] 
      modelYS<-summary(lm(y~matGX[,1]))$coef[2,c(1,2)]
      
      #multiple linear regression (MLR)
      modelXM<-summary(lm(x~matGX[,1]+matC))$coef[2,c(1,2)] #fix for any number of covariates
      modelYM<-summary(lm(y~matGX[,1]+matC))$coef[2,c(1,2)] #fix for any number of covariates
      
      corGXslr<-modelXS[1]/(sqrt( modelXS[1]^2+(n-2)* modelXS[2]^2 ))
      corGXmlr<-modelXM[1]/(sqrt( modelXM[1]^2+(n-2)* modelXM[2]^2 ))
      
      corGYslr<-modelYS[1]/(sqrt( modelYS[1]^2+(n-2)* modelYS[2]^2 ))
      corGYmlr<-modelYM[1]/(sqrt( modelYM[1]^2+(n-2)* modelYM[2]^2 ))
      
      corGX[i]<-corGXmlr - corGXslr
      corGY[i]<-corGYmlr - corGYslr
      
      
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
      
      
      # save for cd-ratio, cd-gls, cd-egger
      u.betaGX_X <- bxU
      u.seGX_X <- sxU
      u.betaGX_Y <- byU
      u.seGX_Y <- syU
      
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
        bxA <- c(bxA, summary(lm(x~matGX[,ll]+matC))$coef[2,1])
        sxA <- c(sxA, summary(lm(x~matGX[,ll]+matC))$coef[2,2])
        
        byA <- c(byA, summary(lm(y~matGX[,ll]+matC))$coef[2,1])
        syA <- c(syA, summary(lm(y~matGX[,ll]+matC))$coef[2,2])
      }
      
      # save for cd-ratio, cd-gls, cd-egger
      a.betaGX_X <- bxA
      a.seGX_X <- sxA
      a.betaGX_Y <- byA
      a.seGX_Y <- syA
      
      for(ll in 1:ncol(matGY)){
        bxA <- c(bxA, summary(lm(x~matGY[,ll]+matC))$coef[2,1])
        sxA <- c(sxA, summary(lm(x~matGY[,ll]+matC))$coef[2,2])
        
        byA <- c(byA, summary(lm(y~matGY[,ll]+matC))$coef[2,1])
        syA <- c(syA, summary(lm(y~matGY[,ll]+matC))$coef[2,2])
      }
      
      
      if("CDcML" %in% colnames.mat1){
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
        
        # CIs - .S.DP
        u.ci.low.XY.S.DP <- res.BDCDcML.u$XtoY.est.S.DP - (qnorm(0.975)*res.BDCDcML.u$XtoY.se.S.DP)
        u.ci.upp.XY.S.DP <- res.BDCDcML.u$XtoY.est.S.DP + (qnorm(0.975)*res.BDCDcML.u$XtoY.se.S.DP)
        
        u.ci.low.YX.S.DP <- res.BDCDcML.u$YtoX.est.S.DP - (qnorm(0.975)*res.BDCDcML.u$YtoX.se.S.DP)
        u.ci.upp.YX.S.DP <- res.BDCDcML.u$YtoX.est.S.DP + (qnorm(0.975)*res.BDCDcML.u$YtoX.se.S.DP)
        
        # decision rules with CIs  --- .S.DP
        if((u.ci.low.XY.S.DP>0 & u.ci.upp.XY.S.DP<1)| (u.ci.low.XY.S.DP> (-1) & u.ci.upp.XY.S.DP<0)){
          if((u.ci.low.YX.S.DP>0 & u.ci.upp.YX.S.DP<1) | (u.ci.low.YX.S.DP> (-1) & u.ci.upp.YX.S.DP<0)){
            u.S.DP <- 3
          }else{
            u.S.DP <- 1
          }
        }else if((u.ci.low.YX.S.DP>0 & u.ci.upp.YX.S.DP<1) & (u.ci.low.YX.S.DP> (-1) & u.ci.upp.YX.S.DP<0)){
          if((u.ci.low.XY.S.DP>0 & u.ci.upp.XY.S.DP<1)| (u.ci.low.XY.S.DP> (-1) & u.ci.upp.XY.S.DP<0)){
            u.S.DP <- 3
          }else{
            u.S.DP <- 2
          }
        }else{
          u.S.DP <- 3
        }
        
        
        
        ############################################
        # Adjusted - CDcML
        ############################################
        a.S.DP<-NA
        
        # CIs - .S.DP
        a.ci.low.XY.S.DP <- res.BDCDcML.a$XtoY.est.S.DP - (qnorm(0.975)*res.BDCDcML.a$XtoY.se.S.DP)
        a.ci.upp.XY.S.DP <- res.BDCDcML.a$XtoY.est.S.DP + (qnorm(0.975)*res.BDCDcML.a$XtoY.se.S.DP)
        
        a.ci.low.YX.S.DP <- res.BDCDcML.a$YtoX.est.S.DP - (qnorm(0.975)*res.BDCDcML.a$YtoX.se.S.DP)
        a.ci.upp.YX.S.DP <- res.BDCDcML.a$YtoX.est.S.DP + (qnorm(0.975)*res.BDCDcML.a$YtoX.se.S.DP)
        
        # decision rules with CIs  --- .S.DP
        if((a.ci.low.XY.S.DP>0 & a.ci.upp.XY.S.DP<1)| (a.ci.low.XY.S.DP>-1 & a.ci.upp.XY.S.DP<0)){
          if((a.ci.low.YX.S.DP>0 & a.ci.upp.YX.S.DP<1) | (a.ci.low.YX.S.DP>-1 & a.ci.upp.YX.S.DP<0)){
            a.S.DP<- 3
          }else{
            a.S.DP<- 1
          }
        }else if((a.ci.low.YX.S.DP>0 & a.ci.upp.YX.S.DP<1) & (a.ci.low.YX.S.DP>-1 & a.ci.upp.YX.S.DP<0)){
          if((a.ci.low.XY.S.DP>0 & a.ci.upp.XY.S.DP<1)| (a.ci.low.XY.S.DP>-1 & a.ci.upp.XY.S.DP<0)){
            a.S.DP<- 3
          }else{
            a.S.DP<- 2
          }
        }else{
          a.S.DP<- 3
        }
        
        
        if(u.S.DP==1 & a.S.DP==1){matR["slr1_mlr1", "CDcML"] <- matR["slr1_mlr1", "CDcML"]+1}
        if(u.S.DP==1 & a.S.DP==2){matR["slr1_mlr2", "CDcML"] <- matR["slr1_mlr2", "CDcML"]+1}
        if(u.S.DP==1 & a.S.DP==3){matR["slr1_mlr3", "CDcML"] <- matR["slr1_mlr3", "CDcML"]+1}
        if(u.S.DP==2 & a.S.DP==1){matR["slr2_mlr1", "CDcML"] <- matR["slr2_mlr1", "CDcML"]+1}
        if(u.S.DP==2 & a.S.DP==2){matR["slr2_mlr2", "CDcML"] <- matR["slr2_mlr2", "CDcML"]+1}
        if(u.S.DP==2 & a.S.DP==3){matR["slr2_mlr3", "CDcML"] <- matR["slr2_mlr3", "CDcML"]+1}
        if(u.S.DP==3 & a.S.DP==1){matR["slr3_mlr1", "CDcML"] <- matR["slr3_mlr1", "CDcML"]+1}
        if(u.S.DP==3 & a.S.DP==2){matR["slr3_mlr2", "CDcML"] <- matR["slr3_mlr2", "CDcML"]+1}
        if(u.S.DP==3 & a.S.DP==3){matR["slr3_mlr3", "CDcML"] <- matR["slr3_mlr3", "CDcML"]+1}
        
      }
      
      
      
      if("MRcML" %in% colnames.mat1){
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
        
        u.p.XY.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.u$XtoY.est.S.DP / res.BDMRcML.u$XtoY.se.S.DP))*2
        u.p.YX.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.u$YtoX.est.S.DP / res.BDMRcML.u$YtoX.se.S.DP))*2
        # decision rules
        if(u.p.XY.BDMRcML.S.DP <= 0.05 & u.p.YX.BDMRcML.S.DP > 0.05){
          u.S.DP<-1
        }#return case 1
        if(u.p.XY.BDMRcML.S.DP > 0.05 & u.p.YX.BDMRcML.S.DP <= 0.05){
          u.S.DP<-2
        }# return case 2
        if((u.p.XY.BDMRcML.S.DP >0.05 & u.p.YX.BDMRcML.S.DP > 0.05) |(u.p.XY.BDMRcML.S.DP <=0.05 & u.p.YX.BDMRcML.S.DP <= 0.05)){
          u.S.DP<-3
        }# return case 3
        
        ############################################
        # Adjusted - MRcCML
        ############################################
        a.S.DP<-NA
        
        a.p.XY.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.a$XtoY.est.S.DP / res.BDMRcML.a$XtoY.se.S.DP))*2
        a.p.YX.BDMRcML.S.DP<- pnorm(-abs(res.BDMRcML.a$YtoX.est.S.DP / res.BDMRcML.a$YtoX.se.S.DP))*2
        # decision rules
        if(a.p.XY.BDMRcML.S.DP <= 0.05 & a.p.YX.BDMRcML.S.DP > 0.05){
          a.S.DP<-1
        }#return case 1
        if(a.p.XY.BDMRcML.S.DP > 0.05 & a.p.YX.BDMRcML.S.DP <= 0.05){
          a.S.DP<-2
        }# return case 2
        if((a.p.XY.BDMRcML.S.DP >0.05 & a.p.YX.BDMRcML.S.DP > 0.05) | (a.p.XY.BDMRcML.S.DP <=0.05 & a.p.YX.BDMRcML.S.DP <= 0.05)){
          a.S.DP<-3
        }# return case 3
        
        if(u.S.DP==1 & a.S.DP==1){matR["slr1_mlr1", "MRcML"] <- matR["slr1_mlr1", "MRcML"]+1}
        if(u.S.DP==1 & a.S.DP==2){matR["slr1_mlr2", "MRcML"] <- matR["slr1_mlr2", "MRcML"]+1}
        if(u.S.DP==1 & a.S.DP==3){matR["slr1_mlr3", "MRcML"] <- matR["slr1_mlr3", "MRcML"]+1}
        if(u.S.DP==2 & a.S.DP==1){matR["slr2_mlr1", "MRcML"] <- matR["slr2_mlr1", "MRcML"]+1}
        if(u.S.DP==2 & a.S.DP==2){matR["slr2_mlr2", "MRcML"] <- matR["slr2_mlr2", "MRcML"]+1}
        if(u.S.DP==2 & a.S.DP==3){matR["slr2_mlr3", "MRcML"] <- matR["slr2_mlr3", "MRcML"]+1}
        if(u.S.DP==3 & a.S.DP==1){matR["slr3_mlr1", "MRcML"] <- matR["slr3_mlr1", "MRcML"]+1}
        if(u.S.DP==3 & a.S.DP==2){matR["slr3_mlr2", "MRcML"] <- matR["slr3_mlr2", "MRcML"]+1}
        if(u.S.DP==3 & a.S.DP==3){matR["slr3_mlr3", "MRcML"] <- matR["slr3_mlr3", "MRcML"]+1}
        
      }
      
      
      
      if("CDRatio" %in% colnames.mat1 | "CDEgger" %in% colnames.mat1 | "CDGLS" %in% colnames.mat1){
        
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
          uR <- 1
        }else if((u.upperCIyx.cdRatio < (-1) | u.lowerCIyx.cdRatio>1) &
                 ((u.lowerCIxy.cdRatio>=(-1) & u.upperCIxy.cdRatio< 0) | (u.lowerCIxy.cdRatio>(0) & u.upperCIxy.cdRatio<= 1)) ){
          uR <- 2
        }else{
          uR <- 3
        }
        
        # CD-Egger decisions
        if((u.upperCIxy.cdEgger < (-1)  | u.lowerCIxy.cdEgger>1) &
           ((u.lowerCIyx.cdEgger>=(-1) & u.upperCIyx.cdEgger<0) | (u.lowerCIyx.cdEgger>(0) & u.upperCIyx.cdEgger<=1))){
          uE <- 1
        }else if((u.upperCIyx.cdEgger < (-1) | u.lowerCIyx.cdEgger>1) &
                 ((u.lowerCIxy.cdEgger>=(-1)  & u.upperCIxy.cdEgger< 0) | (u.lowerCIxy.cdEgger>(0)  & u.upperCIxy.cdEgger<=1))){
          uE <- 2
        }else{
          uE <- 3
        }
        
        # CD-GLS decisions
        if((u.upperCIxy.cdGLS < (-1) | u.lowerCIxy.cdGLS>1) &
           ((u.lowerCIyx.cdGLS>= (-1) & u.upperCIyx.cdGLS<0) | (u.lowerCIyx.cdGLS> (0) & u.upperCIyx.cdGLS<=1)) ){
          uG <- 1
        }else if((u.upperCIyx.cdGLS < (-1) | u.lowerCIyx.cdGLS>1) &
                 ((u.lowerCIxy.cdGLS>= (-1) & u.upperCIxy.cdGLS<0) | (u.lowerCIxy.cdGLS> (0) & u.upperCIxy.cdGLS<=1)) ){
          uG <- 2
        }else{
          uG <- 3
        }
        
        
        ######################################################################################################
        # Adjusted -- CD-Ratio, CD-Egger, CD-GLS (1-Sample, multiple SNPs):
        ######################################################################################################
        
        
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
        
        # CD-RATIO
        if("CDRatio" %in% colnames.mat1){
          aR <- NA
          a.ratio.YX <- a.cd3$CD_Ratio_result$T1toT2
          a.ratio.XY <- a.cd3$CD_Ratio_result$T2toT1
          
          # Confidence intervals for K - needed for use in decision rules
          # CD-ratio CIs
          a.lowerCIyx.cdRatio <- a.ratio.YX["K"] - qnorm(1-sig.level/2)*a.ratio.YX["se(K)"]
          a.upperCIyx.cdRatio <- a.ratio.YX["K"] + qnorm(1-sig.level/2)*a.ratio.YX["se(K)"]
          a.lowerCIxy.cdRatio <- a.ratio.XY["K"] - qnorm(1-sig.level/2)*a.ratio.XY["se(K)"]
          a.upperCIxy.cdRatio <- a.ratio.XY["K"] + qnorm(1-sig.level/2)*a.ratio.XY["se(K)"]
          
          # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
          # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
          # case 3: neither otherwise
          # CD-Ratio decisions
          if(a.upperCIxy.cdRatio < a.lowerCIxy.cdRatio){
            stop("LowerCIxy > UpperCIxy")}
          if(a.upperCIyx.cdRatio < a.lowerCIyx.cdRatio){
            stop("LowerCIyx > UpperCIyx")}
          
          # CD-Ratio decisions
          if((a.upperCIxy.cdRatio < (-1) | a.lowerCIxy.cdRatio>1) &
             ((a.lowerCIyx.cdRatio>=(-1) & a.upperCIyx.cdRatio<0) | (a.lowerCIyx.cdRatio>(0) & a.upperCIyx.cdRatio<=1)) ){
            aR <- 1
          }else if((a.upperCIyx.cdRatio < (-1) | a.lowerCIyx.cdRatio>1) &
                   ((a.lowerCIxy.cdRatio>=(-1) & a.upperCIxy.cdRatio< 0) | (a.lowerCIxy.cdRatio>(0) & a.upperCIxy.cdRatio<= 1)) ){
            aR <- 2
          }else{
            aR <- 3
          }
          
          
          if(uR==1 & aR==1){matR["slr1_mlr1", "CDRatio"] <- matR["slr1_mlr1", "CDRatio"]+1}
          if(uR==1 & aR==2){matR["slr1_mlr2", "CDRatio"] <- matR["slr1_mlr2", "CDRatio"]+1}
          if(uR==1 & aR==3){matR["slr1_mlr3", "CDRatio"] <- matR["slr1_mlr3", "CDRatio"]+1}
          if(uR==2 & aR==1){matR["slr2_mlr1", "CDRatio"] <- matR["slr2_mlr1", "CDRatio"]+1}
          if(uR==2 & aR==2){matR["slr2_mlr2", "CDRatio"] <- matR["slr2_mlr2", "CDRatio"]+1}
          if(uR==2 & aR==3){matR["slr2_mlr3", "CDRatio"] <- matR["slr2_mlr3", "CDRatio"]+1}
          if(uR==3 & aR==1){matR["slr3_mlr1", "CDRatio"] <- matR["slr3_mlr1", "CDRatio"]+1}
          if(uR==3 & aR==2){matR["slr3_mlr2", "CDRatio"] <- matR["slr3_mlr2", "CDRatio"]+1}
          if(uR==3 & aR==3){matR["slr3_mlr3", "CDRatio"] <- matR["slr3_mlr3", "CDRatio"]+1}
          
          
        } # END CD-Ratio section
        
        # START CD-GLS
        if("CDGLS" %in% colnames.mat1){
          aG <- NA
          
          a.gls.YX <- a.cd3$CD_GLS_result$T1toT2
          a.gls.XY <- a.cd3$CD_GLS_result$T2toT1
          
          # CD-GLS CIs - Confidence intervals for K - needed for use in decision rules
          a.lowerCIyx.cdGLS <- a.gls.YX["K"] - qnorm(1-sig.level/2)*a.gls.YX["se(K)"]
          a.upperCIyx.cdGLS <- a.gls.YX["K"] + qnorm(1-sig.level/2)*a.gls.YX["se(K)"]
          a.lowerCIxy.cdGLS <- a.gls.XY["K"] - qnorm(1-sig.level/2)*a.gls.XY["se(K)"]
          a.upperCIxy.cdGLS <- a.gls.XY["K"] + qnorm(1-sig.level/2)*a.gls.XY["se(K)"]
          
          # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
          # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
          # case 3: neither otherwise
          if(a.upperCIxy.cdGLS < a.lowerCIxy.cdGLS){
            stop("LowerCIxy > UpperCIxy")}
          if(a.upperCIyx.cdGLS < a.lowerCIyx.cdGLS){
            stop("LowerCIyx > UpperCIyx")}
          
          # CD-GLS decisions
          if((a.upperCIxy.cdGLS < (-1) | a.lowerCIxy.cdGLS>1) &
             ((a.lowerCIyx.cdGLS>= (-1) & a.upperCIyx.cdGLS<0) | (a.lowerCIyx.cdGLS> (0) & a.upperCIyx.cdGLS<=1)) ){
            aG <- 1
          }else if((a.upperCIyx.cdGLS < (-1) | a.lowerCIyx.cdGLS>1) &
                   ((a.lowerCIxy.cdGLS>= (-1) & a.upperCIxy.cdGLS<0) | (a.lowerCIxy.cdGLS> (0) & a.upperCIxy.cdGLS<=1)) ){
            aG <- 2
          }else{
            aG <- 3
          }
          
          if(uG==1 & aG==1){matR["slr1_mlr1", "CDGLS"] <- matR["slr1_mlr1", "CDGLS"]+1}
          if(uG==1 & aG==2){matR["slr1_mlr2", "CDGLS"] <- matR["slr1_mlr2", "CDGLS"]+1}
          if(uG==1 & aG==3){matR["slr1_mlr3", "CDGLS"] <- matR["slr1_mlr3", "CDGLS"]+1}
          if(uG==2 & aG==1){matR["slr2_mlr1", "CDGLS"] <- matR["slr2_mlr1", "CDGLS"]+1}
          if(uG==2 & aG==2){matR["slr2_mlr2", "CDGLS"] <- matR["slr2_mlr2", "CDGLS"]+1}
          if(uG==2 & aG==3){matR["slr2_mlr3", "CDGLS"] <- matR["slr2_mlr3", "CDGLS"]+1}
          if(uG==3 & aG==1){matR["slr3_mlr1", "CDGLS"] <- matR["slr3_mlr1", "CDGLS"]+1}
          if(uG==3 & aG==2){matR["slr3_mlr2", "CDGLS"] <- matR["slr3_mlr2", "CDGLS"]+1}
          if(uG==3 & aG==3){matR["slr3_mlr3", "CDGLS"] <- matR["slr3_mlr3", "CDGLS"]+1}
        } # end CD-GLS
        
        # START CDEGGER
        if("CDEgger" %in% colnames.mat1){
          aE <- NA
          a.egger.YX<- a.cd3$CD_Egger_result$T1toT2
          a.egger.XY <- a.cd3$CD_Egger_result$T2toT1
          
          # CD-Egger CIs - Confidence intervals for K - needed for use in decision rules
          a.lowerCIyx.cdEgger <- a.egger.YX["K"] - qnorm(1-sig.level/2)*a.egger.YX["se(K)"]
          a.upperCIyx.cdEgger <- a.egger.YX["K"] + qnorm(1-sig.level/2)*a.egger.YX["se(K)"]
          a.lowerCIxy.cdEgger <- a.egger.XY["K"] - qnorm(1-sig.level/2)*a.egger.XY["se(K)"]
          a.upperCIxy.cdEgger <- a.egger.XY["K"] + qnorm(1-sig.level/2)*a.egger.XY["se(K)"]
          
          # case 1: X->Y if CIxy not inside [-1,1] and CIyx inside [-1,0) OR (0,1]
          # case 2: X<-Y if CIyx not inside [-1,1] and CIxy inside [-1,0) OR (0,1]
          # case 3: neither otherwise
          if(a.upperCIxy.cdEgger < a.lowerCIxy.cdEgger){
            stop("LowerCIxyEgger > UpperCIxyEgger")}
          if(a.upperCIyx.cdEgger < a.lowerCIyx.cdEgger){
            stop("LowerCIyxEgger > UpperCIyxEgger")}
          
          # CD-Egger decisions
          if((a.upperCIxy.cdEgger < (-1)  | a.lowerCIxy.cdEgger>1) &
             ((a.lowerCIyx.cdEgger>=(-1) & a.upperCIyx.cdEgger<0) | (a.lowerCIyx.cdEgger>(0) & a.upperCIyx.cdEgger<=1))){
            aE <- 1
          }else if((a.upperCIyx.cdEgger < (-1) | a.lowerCIyx.cdEgger>1) &
                   ((a.lowerCIxy.cdEgger>=(-1)  & a.upperCIxy.cdEgger< 0) | (a.lowerCIxy.cdEgger>(0)  & a.upperCIxy.cdEgger<=1))){
            aE <- 2
          }else{
            aE <- 3
          }
          
          if(uE==1 & aE==1){matR["slr1_mlr1", "CDEgger"] <- matR["slr1_mlr1", "CDEgger"]+1}
          if(uE==1 & aE==2){matR["slr1_mlr2", "CDEgger"] <- matR["slr1_mlr2", "CDEgger"]+1}
          if(uE==1 & aE==3){matR["slr1_mlr3", "CDEgger"] <- matR["slr1_mlr3", "CDEgger"]+1}
          if(uE==2 & aE==1){matR["slr2_mlr1", "CDEgger"] <- matR["slr2_mlr1", "CDEgger"]+1}
          if(uE==2 & aE==2){matR["slr2_mlr2", "CDEgger"] <- matR["slr2_mlr2", "CDEgger"]+1}
          if(uE==2 & aE==3){matR["slr2_mlr3", "CDEgger"] <- matR["slr2_mlr3", "CDEgger"]+1}
          if(uE==3 & aE==1){matR["slr3_mlr1", "CDEgger"] <- matR["slr3_mlr1", "CDEgger"]+1}
          if(uE==3 & aE==2){matR["slr3_mlr2", "CDEgger"] <- matR["slr3_mlr2", "CDEgger"]+1}
          if(uE==3 & aE==3){matR["slr3_mlr3", "CDEgger"] <- matR["slr3_mlr3", "CDEgger"]+1}
        } # End CD-Egger
        
      } # END CD-RATIO / CD-EGGER / CD-GLS 
      
      
    }# end sims loop
    
    
    
    if(plot.pdf){
      
      #########################################
      # Difference in K plot - X --> Y
      #########################################
      cols <- brewer.pal(9, "Set3")
      cols <- cols[c(1:4,9,6:8,5)]
      
      rownames1 <- c("slrC1 mlrC1", "slrC1 mlrC2", "slrC1 mlrC3",
                     "slrC2 mlrC1", "slrC2 mlrC2", "slrC2 mlrC3",
                     "slrC3 mlrC1", "slrC3 mlrC2", "slrC3 mlrC3")
      # barplot
      pdf(paste(plot.name,"_barplot.pdf", sep = ""))
      par(mar=c(5, 4, 4, 6), xpd=TRUE)
      barplot(matR,
              col = cols)
      legend("right", inset = c(-0.2,0), legend = rownames1, 
             fill = cols, box.lty = 0, cex = 0.8,xpd = T)
      dev.off()      

      # boxplot
      pdf(paste(plot.name,"_boxplot.pdf", sep = ""))
      boxplot(corGX, corGY)
      dev.off()     
      
      
      
    } # end plot if statement
    if(table.sv){
      write.table(matR, file=paste0(plot.name,"_seed",SEED,"_matR.txt"), quote=F, row.names=T)
      write.table(corGX, file=paste0(plot.name,"_seed",SEED,"_corGX.txt"), quote=F, row.names=F, col.names = F)
      write.table(corGX, file=paste0(plot.name,"_seed",SEED,"_corGY.txt"), quote=F, row.names=F, col.names = F)
      
    }
    
    
    return(list("matR" = matR, "corGX" = corGX, "corGY"=corGY))
    
    
  }
