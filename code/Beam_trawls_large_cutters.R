##################################################################################
# PULSE TRAWL PROJECT
# Code for Bean large trawls and euro cutters independent scenarios
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

source("~/code/functions.R")

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

## SHRIMP
#csh_beam_small<- cpue_dsvm_input("csh",1)

#setwd("~/Dropbox/PULSE/Pulseproject/data/results")
#files<- mget(ls())
#for (i in 1:length(files)) write.csv(files[[i]], paste(names(files[i]), ".csv", sep = ""))


#Load all .csv files
#----------------------------------------------------------------------------------------
setwd("~/data/results")
temp <- list.files(pattern="*.csv")
for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))

 
sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")

# ---------------------
# ---------------------------------------------
# BEAM TRAWLS: LONG DISTANCE fleet
# --------------------------------------------
# ---------------------
# 
sp1 <- cpue_dsvm_sp(sol_beam_large.csv, sp1, 1.637352)  #SOLE   # gamsol$family$getTheta() for Sole in Jurgens GAMs is 1.637352
sp2 <- cpue_dsvm_sp(ple_beam_large.csv, sp2, 1.195788)  #PLAICE # gample$family$getTheta() for Plaice in Jurgens GAMs is 1.195788
sp3 <- cpue_dsvm_sp(cod_beam_large.csv, sp3, 0.8979249) #COD    # gamcod$family$getTheta() for Cod in Jurgens GAMs is 0.8979249

catchMean(sp4) <- catchMean(sp5) <- array(0.01,dim=dim(catchMean(sp2)),dimnames=dimnames(catchMean(sp2)))
catchSigma(sp4)<- catchSigma(sp5)<- array(0.0000001,dim=dim(catchMean(sp2)),dimnames=dimnames(catchMean(sp2)))

# ---------------------
# PRICES
# ---------------------
price <- read.csv("~/data/input/Visprijzen.csv")

# SOLE = 5 cat 52 weeks3
sp1Price <- price[price$Spec == "Sole" & price$Category != 0 & price$Week != 0,]
sp1Price <- sp1Price[order(sp1Price$Week),]
sp1Price <- array(sp1Price$Price, dim=c(5,52), dimnames=dimnames(catchMean(sp1))[-3])

# Plaice 4 marketcat 1 discard with 0.22 marketvalue (since discarding is not allowed, this fraction goes for fishmeal)
sp2Price <- price[price$Spec == "Plaice" & price$Category != 0 & price$Week != 0,]
# create 5th cat with 0.22 marketvalue (discards)
sp2Price <- rbind(sp2Price,data.frame(Week = 1:52, Category = 5, Price = 0.22, Spec = "Plaice"))
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

effort <- array(c(effort_dsvm_input(4.2,52,11,89.5,"south")[c(1:14,16),8]), dim=c(15,52), dimnames=list(option=dimnames(catchMean(sp1))[[3]],season=as.character(dimnames(catchMean(sp1))[[2]])))

#-------------------------------------------------------------------------------------
# Search for a good sigma value
#-------------------------------------------------------------------------------------

#To find the proper sigma, first I give enough quota to all species and search for a sigma value that allow fishers to go the best 4 or 5 patch options.
#control     <- DynState.control(spp1LndQuota= 160e9,  spp2LndQuota=1000e9, spp1LndQuotaFine= 320, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = 1600, landingCosts= 0.24,gearMaintenance= 87, addNoFishing= TRUE, increments= 20, spp1DiscardSteps= 0, spp2DiscardSteps= 0, sigma= 1.5e+4, simNumber= 1000 , numThreads= 40, verbose=1)

#with sigma equals to    1e+03, probs are 1.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00 by area (1:16 and stay in port)
#with sigma equals to    5e+03, probs are 0.96,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.02,0.00,0.02,0.00,0.00,0.00
#with sigma equals to    1e+04, probs are 0.72,0.00,0.00,0.01,0.00,0.00,0.01,0.00,0.00,0.01,0.10,0.00,0.10,0.02,0.00,0.00 
#with sigma equals to 1.25e+04, probs are 0.61,0.00,0.01,0.02,0.01,0.00,0.02,0.01,0.00,0.02,0.13,0.00,0.13,0.04,0.00,0.00,
#with sigma equals to  1.5e+04, probs are 0.51,0.01,0.02,0.03,0.01,0.01,0.03,0.02,0.00,0.03,0.14,0.00,0.14,0.05,0.00,0.00 
#with sigma equals to    2e+04, probs are 0.38,0.01,0.03,0.05,0.02,0.02,0.04,0.03,0.00,0.05,0.14,0.00,0.15,0.07,0.00,0.00 
#with sigma equals to  2.5e+04, probs are 0.30,0.02,0.04,0.06,0.03,0.03,0.05,0.04,0.00,0.05,0.14,0.01,0.14,0.08,0.00,0.01
#with sigma equals to    1e+05, probs are 0.10,0.05,0.06,0.07,0.06,0.06,0.07,0.06,0.04,0.07,0.09,0.04,0.09,0.07,0.04,0.04


#-------------------------------------------------------------------------------------
# Make contol and execute calculations BEAM TRAWLS with 1600 tons quota
#-------------------------------------------------------------------------------------

control     <- DynState.control(spp1LndQuota= 160e3,  spp2LndQuota=100e9, spp1LndQuotaFine= 320, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = 580, landingCosts= 0.24,gearMaintenance= 87, addNoFishing= TRUE, increments= 20, spp1DiscardSteps= 0, spp2DiscardSteps= 0, sigma= 1.5e+04, simNumber= 1000 , numThreads= 40, verbose=1)



BS160 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_large_B160_south_sigma_1_5e4.RData")

control@spp1LndQuota <- as.numeric(120e3)
BS120 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_large_B120_south_sigma_1_5e4.RData")

control@spp1LndQuota <- as.numeric(80e3)
BS80 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_large_B80_south_sigma_1_5e4.RData")

save(list=c("BS160","BS120","BS80"), file="~/modelresults/Beam_large_south_B160_80_south_sigma_1_5e4.RData.RData")

#-------------------------------------------------------------------------------------
# EFFORT from the NORTH home port; beam trawlers 11 Nm h-1 and 88 mean fishing days                                             
#-------------------------------------------------------------------------------------

effort <- array(c(effort_dsvm_input(6.0,53.4,11,86.4,"north")[c(1:14,16),8]), dim=c(15,52), dimnames=list(option=dimnames(catchMean(sp1))[[3]],season=as.character(dimnames(catchMean(sp1))[[2]])))

# Need to change fuelPrice; fuelconsumption based on relative towing speed while fishing
control     <- DynState.control(spp1LndQuota= 160e3,  spp2LndQuota=100e9, spp1LndQuotaFine= 320, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = 580, landingCosts= 0.24,gearMaintenance= 87, addNoFishing= TRUE, increments= 20, spp1DiscardSteps= 0, spp2DiscardSteps= 0, sigma= 1.5e+04, simNumber= 1000 , numThreads= 40, verbose=1)

BN160 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_large_B160_north_sigma_1_5e4.RData")

control@spp1LndQuota <- as.numeric(120e3)
BN120 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_large_B120_north_sigma_1_5e4.RData")

control@spp1LndQuota <- as.numeric(80e3)
BN80 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_large_B80_north_sigma_1_5e4.RData")

save(list=c("BN160","BN120","BN80"), file="~/modelresults/Beam_large_north_B160_80_sigma_1_5e4.RData")

# ---------------------
# ---------------------------------------------
# SMALL BEAM TRAWLS: EURO CUTTERS fleet
# --------------------------------------------
# here, we need to add the shrimp data
# ---------------------

sp1<- sp2 <- sp3 <- sp4 <- sp5 <-    new("DynStateInput")

# gample$family$getTheta() for Plaice in Jurgens GAMs is 1.195788
# gamsol$family$getTheta() for Sole in Jurgens GAMs is 1.637352
# gamcod$family$getTheta() for Cod in Jurgens GAMs is 0.8979249

sp1 <- cpue_dsvm_sp(sol_beam_small.csv, sp1, 1.637352)  #SOLE
sp2 <- cpue_dsvm_sp(ple_beam_small.csv, sp2, 1.195788)  #PLAICE
sp3 <- cpue_dsvm_sp(cod_beam_small.csv, sp3, 0.8979249) #COD
sp4 <- cpue_dsvm_sp(csh_beam_small.csv, sp4, 1) #SHRIMP in Adriaans GAM family$getTheta() is 1

#although shrimps have the same structure as others, there is only one size class (located in 1)

catchMean(sp5) <- array(0.01,dim=dim(catchMean(sp2)),dimnames=dimnames(catchMean(sp2)))
catchSigma(sp5)<- array(0.0000001,dim=dim(catchMean(sp2)),dimnames=dimnames(catchMean(sp2)))

# ----------------------------------------------------------
# SIZE DEPENDENT PRICING, following Zimmermann et al. (2011)
# -----------------------------------------------------------
# Shrimp 1 marketcat with 1 year marketvalue; strongly affected by the volume landed
# we pick up the mean price from 2010-2015
# add some variability dependant on the observed cpue

price<- mean(subset(read.csv("~/data/input/priceshrimp.csv"), year %in% c(2010:2015))$value)
slopeprice<- 0.01
wts<- csh_beam_small.csv$data
wts[wts<0]<-0

pred      <- csh_beam_small.csv
pred$data <- c(price + slopeprice*(((wts-mean(wts))/mean(wts))))
pred <- pred[,c("sizeclass","data","area", "season")]
names (pred) <- c("cat", "data", "option", "season" )
pred <- tapply(pred[,"data"], list(cat=factor(x = pred[,"cat"], levels = unique(pred[,"cat"])), 
                                     season=factor(x = pred[,"season"], levels = as.character(sort(unique(pred[,"season"])))), 
                                     option=factor(x = pred[,"option"], levels = unique(pred[,"option"]))), sum)
pred[2:5,,]<-0  # remove prices for other categories different than 1
#select the area 2 where there is more variability in catches and therefore in prices
sp4Price <- array(c(pred[,,2]), dim=c(5,52), dimnames=dimnames(catchMean(sp3))[-3])

#-------------------------------------------------------------------------------------
# EFFORT from the SOUTH home port; euro cutters 8.9 Nm h-1 and 88 mean fishing days                                               
#-------------------------------------------------------------------------------------

effort <- array(c(effort_dsvm_input(4.2,52,8.9,89.5,"south")[c(1:14,16),8]), dim=c(15,52), dimnames=list(option=dimnames(catchMean(sp1))[[3]],season=as.character(dimnames(catchMean(sp1))[[2]])))

# Need to change fuelPrice; fuelconsumption based on relative towing speed while fishing
control     <- DynState.control(spp1LndQuota= 160e3,  spp2LndQuota=100e9, spp1LndQuotaFine= 320, spp2LndQuotaFine= 3e6, fuelUse = 1, fuelPrice = 162, landingCosts= 0.24,gearMaintenance= 87, addNoFishing= TRUE, increments= 20, spp1DiscardSteps= 0, spp2DiscardSteps= 0, sigma= 1.5e+04, simNumber= 1000 , numThreads= 40, verbose=1)
 
CS160 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_small_C160_south_sigma_1_5e4.RData")

control@spp1LndQuota <- as.numeric(120e3)
CS120 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_small_C120_south_sigma_1_5e4.RData")

control@spp1LndQuota <- as.numeric(80e3)
CS80 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_small_C80_south_sigma_1_5e4.RData")



save(list=c("CS160","CS120","CS80"), file="~/modelresults/Beam_small_south_C160_80_sigma_1_5e4.RData")

#-------------------------------------------------------------------------------------
# EFFORT from the NORTH home port; euro cutters 8.9 Nm h-1 and 88 mean fishing days                                                      
#-------------------------------------------------------------------------------------

effort <- array(c(effort_dsvm_input(6.0,53.4,8.9,86.4,"north")[c(1:14,16),8]), dim=c(15,52), dimnames=list(option=dimnames(catchMean(sp1))[[3]],season=as.character(dimnames(catchMean(sp1))[[2]])))

control@spp1LndQuota <- as.numeric(160e3)
CN160 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_small_B160_north_sigma_1_5e4.RData")

control@spp1LndQuota <- as.numeric(120e3)
CN120 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_small_B120_north_sigma_1_5e4.RData")

control@spp1LndQuota <- as.numeric(80e3)
CN80 <- DynState(sp1, sp2, sp3, sp4, sp5, sp1Price, sp2Price, sp3Price, sp4Price, sp5Price, effort, control)
save.image("~/modelresults/Beam_small_B80_north.RData")

save(list=c("CN160","CN120","CN80"), file="~/modelresults/Beam_small_north_C160_80_sigma_1_5e4.RData")






# control <- DynState.control(Increments=30, PlaiceUplimit=1600000,PlaiceDiscardSteps=1, SoleDiscardSteps=0,CodDiscardSteps=0, 
#                             EffortUplimit=NA,Handling= 0.24,CrewShare = 0.33,GearCost=87,VarCost= 0.05,EffortPrice=600,
#                             FinePlaice=320,ChoiceDist=1,SimNumber=1000, NumThreads=20)

#area<- c(10:15)
#time<- 1:52
#sp1@catchMean <- sp1@catchMean[,time,area]
#sp1@catchSigma<- sp1@catchSigma[,time,area]
#sp2@catchMean <- sp2@catchMean[,time,area]
#sp2@catchSigma<- sp2@catchSigma[,time,area]
#sp3@catchMean <- sp3@catchMean[,time,area]
#sp3@catchSigma<- sp3@catchSigma[,time,area]
#sp4@catchMean <- sp4@catchMean[,time,area]
#sp4@catchSigma<- sp4@catchSigma[,time,area]
#sp5@catchMean <- sp5@catchMean[,time,area]
#sp5@catchSigma<- sp5@catchSigma[,time,area]

#sp1Price<- sp1Price[,time]
#sp2Price<- sp2Price[,time]
#sp3Price<- sp3Price[,time]
#sp4Price<- sp4Price[,time]
#sp5Price<- sp5Price[,time]

#effort<- effort[area, time]160e3 160000

