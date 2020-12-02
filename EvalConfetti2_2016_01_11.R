# Compiling under windows
#R CMD SHLIB -o Confetti2_01_12_norm.dll R_Main.cpp Forest.cpp mersenne.cpp stoc1.cpp  userintf.cpp

# R --arch x64 CMD SHLIB -o Confetti2_01_11_norm.dll R_Main.cpp Forest.cpp mersenne.cpp stoc1.cpp  userintf.cpp

# Compiling under Linux
# module load R/3.1.0-1_gcc_4.8.1
# module load RStudio/0.98.507_R-3.1.0
# R CMD SHLIB -o Confetti2_2015_10_30.so R_Main.cpp Forest.cpp mersenne.cpp stoc1.cpp  userintf.cpp

#---------------------------------------------------------------------------------
Confetti2 <- function(parvec,
                      nGen=100,
                      nRep=1,
                      forest="BCI"
)                
{
   # BCI settings --------------------------------------------------
   if (forest == "BCI"){
      xmax <- 1000
      ymax <- 500
      nTrees <- 21000
      Jmeta <- 2e6
      map.cell.size <- 20 # [m]
      
      map.file <- "InOut\\HabitatMapBCI_Harms_nomixed.txt"
      rel.dens.file <- "InOut\\RelativeDensityBCI_harms_nomixed_n50.txt"
      
      #map.file <- "InOut/HabitatMapBCI_Harms_nomixed.txt"
      #rel.dens.file <- "InOut/RelativeDensityBCI_harms_nomixed_n50.txt"
      
      n.hab.types <- 6
   }
   # ---------------------------------------------------------------
   else {
   
      # Sinharaja settings --------------------------------------------
      if (forest == "Sin"){
         xmax <- 500
         ymax <- 500
         nTrees <- 17500
         Jmeta <- 2e6
         map.cell.size <- 500/26
         
         map.file <- "InOut\\HabitatMapSinharaja_Ruwan1.txt"
         rel.dens.file <- "InOut\\RelativeDensitySinharaja_Ruwan_n50.txt"
         
         #map.file <- "InOut/HabitatMapSinharaja_Ruwan1.txt"
         #rel.dens.file <- "InOut/RelativeDensitySinharaja_Ruwan_n50.txt"
         
         n.hab.types <- 5
      }  
      # ---------------------------------------------------------------
      else return ("Error in forest label!")
   }   
   
   #so.file <- "Confetti2_12_07_norm.dll"
   so.file <- "Confetti2_01_11_norm.dll"
   
   dyn.load(so.file)

   #seed1 <- parvec[1]
   
   parvec <- as.numeric(parvec)
   
   #parvec[c(1,2,4,6:9)] <- 10^parvec[c(1,2,4,6:9)]
   #print(paste("\nParameters:",paste(parvec[c(1,2,7)],collapse=" ","\n")))
   #print(parvec)
   
   NClassSAD <- 12
   
   if (xmax>ymax) NScalesSAR <- 13
   else           NScalesSAR <- 12 
   
   nspec.max <- 1e4
   nspecpairs.max <- 2e5 
   
   abund.min <- 50
   rmax1     <- 50
   
   # output data
   Seeds     <- numeric(nRep)
   PredNSpec <- numeric(nRep)
   PredShannon <- numeric(nRep)
   PredAnnualMort <- numeric(nRep)
   PredBD_5 <- numeric(nRep)
   PredBD_tot <- numeric(nRep)
   
   PredSAD <- matrix(0,nrow=nRep,ncol=NClassSAD)
   PredPCF <- matrix(0,nrow=nRep,ncol=rmax1)
   PredPropCon <- matrix(0,nrow=nRep,ncol=rmax1)
   PredSAR1 <- matrix(0,nrow=nRep,ncol=rmax1)
   PredSAR2_m <- matrix(0,nrow=nRep,ncol=NScalesSAR)
   PredSAR2_sd <- matrix(0,nrow=nRep,ncol=NScalesSAR)
   
   Pred.mMi <- matrix(NA,nrow=nRep,ncol=2)
   Pred.sdMi <- matrix(NA,nrow=nRep,ncol=2)
   colnames(Pred.mMi) <- colnames(Pred.sdMi) <- c("10m","50m")
   
   Pred.mMij <- matrix(NA,nrow=nRep,ncol=2)
   Pred.sdMij <- matrix(NA,nrow=nRep,ncol=2)
   colnames(Pred.mMij) <- colnames(Pred.sdMij) <- c("10m","50m")
      
   for (irep in 1:nRep){

      model.pred <- .C("PredFDP",
                      nGen = as.integer(nGen),
                      #seed = seed1,
                      seed = as.integer(99),
                      #seed = as.integer(runif(1,1,999999)),
                      xmax = as.double(xmax),
                      ymax = as.double(ymax),
                      ntrees = as.integer(nTrees),
                      jmeta = as.integer(Jmeta),
                      map_cell_size = as.double(map.cell.size),
                      map_file = as.character(map.file),
                      rel_dens_file = as.character(rel.dens.file),
                      n_hab_types = as.integer(n.hab.types),
                      # model parameters
                      theta = as.double(parvec[1]),
                      m     = as.double(parvec[2]),
                      Rmax = as.double(parvec[3]),
                      aRec = as.double(parvec[4]),
                      aHab = as.double(parvec[5]),
                      aSurv = as.double(parvec[6]),
                      bSurv = as.double(parvec[7]),
                      m_dm_spec = as.double(parvec[8]),
                      sd_dm_spec = as.double(parvec[9]),
                      m_JC = as.double(parvec[10]),
                      sd_JC = as.double(parvec[11]),
                      #output
                      NSpec = as.integer(0),
                      Shannon = as.double(0),
                      Abund = as.integer(numeric(nspec.max)),
                      NClass = as.integer(NClassSAD),
                      SAD = as.integer(numeric(NClassSAD)),
                      AnnualMort = as.double(0),
                      BD_5 = as.integer(0),
                      BD_tot = as.integer(0),
                      PCF = as.double(numeric(rmax1)),
                      PropCon = as.double(numeric(rmax1)),
                      SAR1 = as.double(numeric(rmax1)),
                      SAR2_m = as.double(numeric(NScalesSAR)),
                      SAR2_sd = as.double(numeric(NScalesSAR)),
                      # species specific patterns
                      Kcon10 = as.double(numeric(nspec.max)),
                      Kcon50 = as.double(numeric(nspec.max)),
                      Khet10 = as.double(numeric(nspecpairs.max)),
                      Khet50 = as.double(numeric(nspecpairs.max)) #,
                      #Dhet10 = as.double(numeric(nspecpairs.max)),
                      #Dhet50 = as.double(numeric(nspecpairs.max)),
                      #AbundSpec2 = as.integer(numeric(nspecpairs.max)),
                      #xPOD = as.double(numeric(nspecpairs.max)) ,
                      #NNDist = as.double(numeric(nspecpairs.max)),
                      # tree specific data
                      #gx = as.double(numeric(nTrees)),
                      #gy = as.double(numeric(nTrees)),
                      #sp = as.integer(numeric(nTrees))            
      )
      
      Seeds[irep]     <- model.pred$seed
      
      PredNSpec[irep] <- model.pred$NSpec
      PredShannon[irep] <- model.pred$Shannon
      PredBD_5[irep] <- model.pred$BD_5
      PredBD_tot[irep] <- model.pred$BD_tot
      PredAnnualMort[irep] <- model.pred$AnnualMort

      PredSAD[irep,] <- model.pred$SAD
      PredPCF[irep,] <- model.pred$PCF
      PredPropCon[irep,] <- model.pred$PropCon
      PredSAR1[irep,] <- model.pred$SAR1
    
      PredSAR2_m[irep,] <- model.pred$SAR2_m
      PredSAR2_sd[irep,] <- model.pred$SAR2_sd
      
      abund <- model.pred$Abund[model.pred$Abund>0]
      abund2 <- model.pred$Abund[model.pred$Abund>=abund.min]
            
#       logKcon10 <- log10(model.pred$Kcon10[1:length(abund2)])
#       logKcon50 <- log10(model.pred$Kcon50[1:length(abund2)])
#        
#       logKhet10 <- log10(model.pred$Khet10[model.pred$Khet10>0])
#       logKhet50 <- log10(model.pred$Khet50[model.pred$Khet50>0])

      Mi10 <- log10(model.pred$Kcon10[1:length(abund2)]) - log10(pi*10^2)
      Mi50 <- log10(model.pred$Kcon50[1:length(abund2)]) - log10(pi*50^2)
             
      Mij10 <- log10(model.pred$Khet10[model.pred$Khet10>0]) - log10(pi*10^2)
      Mij50 <- log10(model.pred$Khet50[model.pred$Khet50>0]) - log10(pi*50^2)
 
#       nspec.pairs <- length(abund2)^2 - length(abund2)
#       Dhet10 <- model.pred$Dhet10[1:nspec.pairs]
#       Dhet50 <- model.pred$Dhet50[1:nspec.pairs]
# 
#       bivariate emptiness probability
#       see Wiegand et al. 2007 AmNat
#       dens2 <- model.pred$AbundSpec2[1:nspec.pairs] /(1000*500) 
#       Phet10 <- -(1-Dhet10) + exp(-dens2*pi*10^2)
#       Phet50 <- -(1-Dhet50) + exp(-dens2*pi*50^2)
# 
#       Chet10 <- -log(1-Dhet10)/(dens2*pi*10^2)
#       Chet50 <- -log(1-Dhet50)/(dens2*pi*50^2)
#       
#       xPOD <- model.pred$xPOD[abs(model.pred$xPOD)>0]
#        
#        NNDist <- model.pred$NNDist[model.pred$NNDist>0]
#        NNDist <- NNDist[NNDist<1500.0]
#       
#       #Negative Binomial ditribution
#       cens1 <- data.frame(gx=model.pred$gx,
#                           gy=model.pred$gy,
#                           sp=as.factor(model.pred$sp))
#       clump1 <- GetSpecNBD(cens1,plot.size=50)
#       log.kNB <- log10(clump1[,"kNB"])
#       Psq0 <- clump1[,"p0"]
#      
      Pred.mMi[irep,1]  <- mean(Mi10[is.finite(Mi10)])
      Pred.sdMi[irep,1] <- sd(Mi10[is.finite(Mi10)])
      Pred.mMi[irep,2]  <- mean(Mi50[is.finite(Mi50)])
      Pred.sdMi[irep,2] <- sd(Mi50[is.finite(Mi50)])

      Pred.mMij[irep,1]  <- mean(Mij10[is.finite(Mij10)])
      Pred.sdMij[irep,1] <- sd(Mij10[is.finite(Mij10)])
      Pred.mMij[irep,2]  <- mean(Mij50[is.finite(Mij50)])
      Pred.sdMij[irep,2] <- sd(Mij50[is.finite(Mij50)])
        
#       Pred.mxPOD[irep] <- mean(xPOD)
#       Pred.sdxPOD[irep] <- sd(xPOD) 
#       
#       Pred.mkNB[irep] <- mean(log.kNB)
#       Pred.sdkNB[irep] <- sd(log.kNB)
#       
#       Pred.mPsq0 <- mean(Psq0)
#       Pred.sdPsq0 <- sd(Psq0)
   }
    
   dyn.unload(so.file)

   out1 <- list()
   
   out1$Seeds <- Seeds
   out1$pars <- parvec
   
   out1$mNSpec <- mean(PredNSpec)
   out1$mShannon <- mean(PredShannon)
   out1$mBD_5 <- mean(PredBD_5)
   out1$mBD_tot <- mean(PredBD_tot)
   out1$mAnnualMort <- mean(PredAnnualMort)
   
   out1$mSAD <- colMeans(PredSAD)

   out1$mPCF <-  colMeans(PredPCF)
   out1$mPropCon <- colMeans(PredPropCon)
   out1$mSAR1 <- colMeans(PredSAR1)
   out1$mSAR2_m <- colMeans(PredSAR2_m)
   out1$mSAR2_sd <- colMeans(PredSAR2_sd)
   
   out1$mMi <- colMeans(Pred.mMi)
   out1$sdMi <- colMeans(Pred.sdMi)
    
   out1$mMij <- colMeans(Pred.mMij)
   out1$sdMij <- colMeans(Pred.sdMij)
    
#   out1$mxPOD <- mean(Pred.mxPOD)
#   out1$sdxPOD <- mean(Pred.sdxPOD)
   
#   out1$mkNB <- mean(Pred.mkNB)
#   out1$sdkNB <- mean(Pred.sdkNB)
    
#   out1$mPsq0 <- mean(Pred.mPsq0)
#   out1$sdPsq0 <- mean(Pred.sdPsq0)
   
   out1$Mi10 <- Mi10
   out1$Mi50 <- Mi50
    
   out1$Mij10 <- Mij10
   out1$Mij50 <- Mij50
 
#   out1$Dhet10 <- Dhet10
#   out1$Dhet50 <- Dhet50
# 
#   out1$Chet10 <- Chet10
#   out1$Chet50 <- Chet50
# 
#   out1$xPOD    <- xPOD
    out1$Abund <- abund
    out1$Abund2 <- abund2

#    out1$NNDist <- NNDist
#    out1$logkNB <- log.kNB
#    out1$Psq0 <- Psq0
   
    return(out1)
}


#Test Confetti
#read data --------------------------------------------------------------------
#   
para1 <- c(theta=50,
           m=0.1,
           Rmax=5.0,
           aRec=1.0,
           aHab=1.0,
           aSurv=999,
           bSurv=0.89,
           m_dm_spec=30,
           sd_dm_spec=0,
           m_JCspec=1,
           sd_JC_spec=0
           )

#call function
system.time(check <- Confetti2(parvec = para1,nRep=1,nGen=100,forest="BCI"))
str(check)
check$mBD_tot/21000
plot(check$mPropCon,type="l",ylim=c(0,0.1),lwd=2)

propcon2 <- as.numeric(read.table("InOut/PropCon1.csv",header=T,sep=";"))
points(1:50,propcon2,col="red",lty=2,lwd=2)
