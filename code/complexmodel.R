##################################################################################
# PULSE TRAWL PROJECT
# 8 October 2018
# IJmuiden- IMARES
##################################################################################

options(width=240)
library(RDynState5NAsigmaseason52Age5)
library(ggplot2)
library(reshape2)
library(FLCore)
library(plyr)
library(scales)

source("~/Dropbox/PULSE/Pulseproject/code/functions.R")

#----------------------------------------------------------------------------------------
# Rerun GAM prediction with updated catch efficiency parameters
#----------------------------------------------------------------------------------------

## PLAICE
#ple_beam_large<- cpue_dsvm_input("ple",1)
#ple_pulse_large<- cpue_dsvm_input("ple",1.28)
#ple_beam_small<- cpue_dsvm_input("ple",(1*0.282))
#ple_pulse_small<- cpue_dsvm_input("ple",(1.28*0.243))

## SOLE
#sol_beam_large<- cpue_dsvm_input("sol",1)
#sol_pulse_large<- cpue_dsvm_input("sol",0.6)
#sol_beam_small<- cpue_dsvm_input("sol",(1*0.282))
#sol_pulse_small<- cpue_dsvm_input("sol",(0.6*0.243))

## COD
#cod_beam_large<- cpue_dsvm_input("cod",1)
#cod_pulse_large<- cpue_dsvm_input("cod",0.6)
#cod_beam_small<- cpue_dsvm_input("cod",(1*0.282))
#cod_pulse_small<- cpue_dsvm_input("cod",(0.6*0.243))

#setwd("~/Dropbox/PULSE/Pulseproject/data/results")
#files<- mget(ls())
#for (i in 1:length(files)) write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))


#Load all .csv files
#----------------------------------------------------------------------------------------
setwd("~/Dropbox/PULSE/Pulseproject/data/results")
temp <- list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

 
sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")

# gample$family$getTheta() for Plaice in Jurgens GAMs is 1.195788
# gamsol$family$getTheta() for Sole in Jurgens GAMs is 1.637352
# gamcod$family$getTheta() for Cod in Jurgens GAMs is 0.8979249

# ---------------------
# Bean trawls long distance fleet
# ---------------------

sp1 <- cpue_dsvm_sp(sol_beam_large.csv, sp1, 1.637352)  #SOLE
sp2 <- cpue_dsvm_sp(ple_beam_large.csv, sp2, 1.195788)  #PLAICE
sp3 <- cpue_dsvm_sp(cod_beam_large.csv, sp3, 0.8979249) #COD

catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=dim(catchMean(sp2)),dimnames=dimnames(catchMean(sp2)))
catchSigma(sp4)<- catchSigma(sp5)<- array(0.0000001,dim=dim(catchMean(sp2)),dimnames=dimnames(catchMean(sp2)))

# ---------------------
# SIZE DEPENDENT PRICING, following Zimmermann et al. (2011)
# ---------------------
price <- read.csv("~/Dropbox/PULSE/Pulseproject/data/input/Visprijzen.csv")

# sp1Price <- array(c(sp1price + slope1price*(((wts-mean(wts))/mean(wts)))), dim=c(length(ages),length(season)), dimnames=list(cat=ages,season=as.character(season)))

# SOLE = 5 cat 52 weeks3
sp1Price <- price[price$Spec == "Sole" & price$Category != 0 & price$Week != 0,]
sp1Price <- sp1Price[order(sp1Price$Week),]
sp1Price <- array(sp1Price$Price, dim=c(5,52), dimnames=dimnames(catchMean(sp1))[-3])

# Plaice 4 marketcat 1 discard with 0 marketvalue
sp2Price <- price[price$Spec == "Plaice" & price$Category != 0 & price$Week != 0,]
# create 5th cat with no marketvalue (discards)
sp2Price <- rbind(sp2Price,data.frame(Week = 1:52, Category = 5, Price = 0, Spec = "Plaice"))
sp2Price <- sp2Price[order(sp2Price$Week),]
sp2Price <- array(sp2Price$Price, dim=c(5,52), dimnames=dimnames(catchMean(sp2))[-3])

# Cod 5 marketcat with 52 weeks marketvalue
sp3Price <- price[price$Spec == "Cod" & price$Category != 0 & price$Category != 6 & price$Week != 0,]      #delete cat 0 and 6 and no week 0
sp3Price <- sp3Price[order(sp3Price$Week),]
sp3Price <- array(sp3Price$Price, dim=c(5,52), dimnames=dimnames(catchMean(sp3))[-3])

sp4Price <- sp5Price <- array(c(0), dim=c(5,52), dimnames=dimnames(catchMean(sp1))[-3])

#-------------------------------------------------------------------------------------
# EFFORT from the SOUTH home port; beam trawlers 11 Nm h-1 and 88 mean fishing days                                                      
#-------------------------------------------------------------------------------------
#URK
#effort  <- array(c(14,16,14,14,15,17,15,15,18),dim=9,2)       # kwart van een dag, so total effort divided by 4!

effort <- array(c(effort_dsvm_input(4.2,52,11,89.5,"south")[c(1:14,16),8]), dim=c(15,52), dimnames=list(option=dimnames(catchMean(sp1))[[3]],season=as.character(dimnames(catchMean(sp1))[[2]])))

#-------------------------------------------------------------------------------------
# Make contol and execute calculations BEAM TRAWLS with 1600 tons quota
#-------------------------------------------------------------------------------------
# control <- DynState.control(Increments=30, PlaiceUplimit=1600000,PlaiceDiscardSteps=1, SoleDiscardSteps=0,CodDiscardSteps=0, 
#                             EffortUplimit=NA,Handling= 0.24,CrewShare = 0.33,GearCost=87,VarCost= 0.05,EffortPrice=600,
#                             FinePlaice=320,ChoiceDist=1,SimNumber=1000, NumThreads=20)
control     <- DynState.control(spp1LndQuota= 160000,  spp2LndQuota=1000000, spp1LndQuotaFine= 320, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = 1600, landingCosts= 0,gearMaintenance= 0, addNoFishing= TRUE, increments= 20, spp1DiscardSteps= 0, spp2DiscardSteps= 0, sigma= 1, simNumber= 1000 , numThreads= 20, verbose=1)
 

z <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)





effort_dsvm_input(6.0,53.4,11,88,"north")
