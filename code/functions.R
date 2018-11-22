##################################################################################################
# FUNCTIONS FOR THE PULSE PROJECT
# 23rd October 2018
# Nekane Alzorriz
###################################################################################################

library(mgcv)
library(lattice)
library(fields)

#---------------------------------------------------------------------------------------------------------
# Effort calculations for model-runs
#---------------------------------------------------------------------------------------------------------


effort_dsvm_input <- function(lon, lat, steamspeed, fishingtime, port){
 
#---------------------------------------------------------------------------------------------------------
# determine mean fishing-effort large beam trawlers or euro-cutters
#---------------------------------------------------------------------------------------------------------
# need the lon and lat from the home port to estimate the distance
# need the mean steam speed and fishing time of the studied fleet

# from home port
dist        <- read.csv("~/data/input/coord_distance.csv")
dist        <- dist[dist$port==as.character(port),]
dist        <- dist[,c(3,2)]
area        <- dist[1:16,]
homeport    <- dist[16,]
homeport[]  <- c(lon, lat)

dist_to_port <- rdist.earth(area, homeport, miles = F, R = NULL)*2                    # times 2 for return travel!
colnames(dist_to_port)[1] <- "km_dist"

# conversions-factor 1Nm = 1.825km
Nm           <- rep(1.852,16)                                                            
dist_to_port <- cbind(dist_to_port,Nm)
Nm_dist      <- dist_to_port[,1]/dist_to_port[,2]
newdat       <- cbind(dist_to_port,Nm_dist)

# vessels have steamingspeed of 12Nm h-1 (nav stuk mike en Jan Jaap)
speed  <- rep(steamspeed,16)                                                             
newdat <- cbind(newdat,speed)

# Calculate steam-time to area
steamT              <- newdat[,3]/speed
newdat              <- cbind(newdat,steamT)
colnames(newdat)[5] <- "steam_time"

# Calculate tot time per trip = steam + fishing time (hours)
fishing  <- rep(fishingtime,16)
tot_time <- fishing + newdat[,5]                                               
newdat   <- cbind(newdat,fishing,tot_time)  

# calculate Days-at-Sea 
#DAS1<- newdat[,7]/24
#DAS <- round(newdat[,7]/12)
#newdat <- cbind(newdat,DAS)

# 3 2 3 3 2 2 3 3 3 2 3 3 3 3 4 3 Low variation among areas!!!

# to get more differentiation between the areas 
# devided tot.time by 6 (quarter of a day)
DAS_UK <- ceiling(newdat[,7]/6)
newdat <- cbind(newdat,DAS_UK)

return(newdat)

}


#---------------------------------------------------------------------------------------------------------
# Reestimate catches for model runs based on catch efficiency and stock biomass corrections
#---------------------------------------------------------------------------------------------------------


cpue_dsvm_input <- function(sp, catchefficiency){
  
  #---------------------------------------------------------------------------------------------------------
  # Load GAM data
  #---------------------------------------------------------------------------------------------------------
  #  RUN WITH R2.14.2 and mgcv 1.7.13
  # need the lon and lat from the home port to estimate the distance
  # need the mean steam speed and fishing time of the studied fleet
  
  # COD
  if (sp=="cod"){
    load("~/data/input/GAMdata/gamcod_final.rdata")
    gamdata <- gamcod
    # Use Age-Length_key mean weight per month (data-file for points where no month-area match is observed)
    Key     <- read.csv("~/Dropbox/PULSE/Pulseproject/data/input/snij_Codb.csv")
    stockcorrection <- 4.38 # cod conversion factor from 70s to 2000
  }
  # PLAICE
  if (sp=="ple"){
    load("~/data/input/GAMdata/gample_final.rdata")
    gamdata <- gample
    Key     <- read.csv("~/data/input/snij_Pleb.csv")
    disc    <- read.csv("~/data/input/discards_in_weight.csv")
    disc$X  <- NULL
    colnames(disc)[3]<- "data"
    disc    <- disc[disc$week!= 53,] # remove plaice box areas + area 1...JUST REMOVE WEEK 53
    stockcorrection <- 0.78 #plaice, but only for what is not discards, !=5 sizeclass
  }
  # SOLE
  if (sp=="sol"){
    load("~/data/input/GAMdata/gamsol_final.rdata")
    gamdata <- gamsol
    Key     <- read.csv("~/data/input/snij_Solb.csv")
    stockcorrection <- 1.23 # sole conversion factor from 70s to 2000
  }
  
  if(sp!="csh"){
    # use expand.grid() because of multiple size classes and weeks (16 areas * 5 classes * 52 weeks = 4160)
    pred <- expand.grid(sizeclass=1:5,area=as.factor(c(1:16)),week=1:52)
    pred <- cbind(pred, gear="BT", lhp=log(2000), year=1974, yrwk=1, fishhour = 88)
    pred$yearwk <- (pred$year - 1970)*52 + pred$week             # vormt extra kolom!
    pred$data<-rep(NA,nrow(pred))
    
    # Do prediction
    pred_with_se   <- predict(gamdata,newdata=pred,se=T) #get predicted value
    pred$data   <- exp(pred_with_se$fit)
    pred$season <- pred$week
    pred$month  <- ceiling(pred$week/4.34)
    
    # make area 12 equal to area 16 
    pred[pred$area == 12,]$data <- pred[pred$area == 16,]$data
    
    # Use Age-Length_key
    dat <- expand.grid(month = 1:12, area=c(1:16), sizeclass=1:5)
    dat <- merge(dat,Key, by=c("month","sizeclass"), all.x=T)
    dat <- dat[,c(1,2,3,4)]
    
    # merge predictions with mean weight per month and area
    pred  <- merge(pred,dat, by=c("sizeclass","area","month"))            # merge predictions with weight-file
    pred  <- pred[order(pred$area),]
    pred$weight_kg <- pred$mean_weight/1000                                   # zet om van gr -> kg
    pred$data.Weight <- pred$data * pred$weight_kg                         # calculate data in weight(kg) per trip
    
    # select necessary columns
    pred        <- pred[,c(1,15,2,12)]                                             # select columns to fit pred.input array (sizeclass, data, area, season)
    #pred$data   <- pred$data.Weight - pred$data.Weight * (1- catchefficiency)      # converge BT or pulse data based on catch efficiency values (69.3% kg/h less cod paper Bob v. Marlen)
    pred$data   <- pred$data.Weight * catchefficiency 
    pred        <- pred[,c(1,5,3,4)]                                               # select columns to fit input array (sizeclass, data, area, season)
    pred        <- pred[order(pred$season),]                                       # order to get right sequence in pred.input
    if (sp=="ple") pred[is.na(pred$data),]$data <- disc$data                       # enter discard data
    
    # combine area 15 and 16 by averaging predictions
    area16  <- pred[pred$area %in% c(15,16),]
    area16  <- aggregate(area16$data, by=list(sizeclass=area16$sizeclass, season=area16$season),mean)
    colnames(area16)[3]<-"data"
    area16$area<- 16
    area16  <- area16[,c(1,3,4,2)]
    
    # deselect area 15 & 16 and fill with new results area
    pred <- pred[!pred$area %in% c(15,16),]
    pred<- rbind(pred,area16)
    # change data to resemble current stocks (SSB 2008-2012)
    #  if (sp =="ple"){
    #       pred[pred$sizeclass != 5,]$data <- pred[pred$sizeclass != 5,]$data/stockcorrection
    #      } else {
    #        pred$data <- pred$data/stockcorrection
    #      }
    pred$data <- pred$data/stockcorrection
  }
  
  if(sp=="csh"){
    gamdata <- read.csv("~/data/input/GAMdata/cpue_csh.csv")
    gamdata$sizeclass<- 1
    names(gamdata)[2]<- "season"
    names(gamdata)[5]<- "data"
    gamdata <- gamdata[,c(1,8,5,3,2)] 
    dat <- expand.grid(sizeclass=1:5,area=as.factor(c(1:16)),season=1:52)
    pred  <- merge(dat,gamdata,by=c("sizeclass","area","season"),all.x = TRUE)            # merge predictions with weight-file
    pred  <- pred[order(pred$area),]
    pred[is.na(pred)]<-0
    pred <- pred[!pred$area %in% c(15),]
    # select necessary columns
    pred <- pred[,c(4,1,5,2,3)]                                             # select columns to fit pred.input array (sizeclass, data, area, season)
    
  }
    
  return(pred)
  
}

#---------------------------------------------------------------------------------------------------------
# GAM cpue data to data frame and specified in wich sp DSVM location want to allocate
#---------------------------------------------------------------------------------------------------------

cpue_dsvm_sp <- function(data, sp, theta){

    pred <- data[,c("sizeclass","data","area", "season")]
    names (pred) <- c("cat", "data", "option", "season" )
    catchMean(sp)  <- tapply(pred[,"data"], list(cat=factor(x = pred[,"cat"], levels = unique(pred[,"cat"])), 
                                                   season=factor(x = pred[,"season"], levels = as.character(sort(unique(pred[,"season"])))), 
                                                   option=factor(x = pred[,"option"], levels = unique(pred[,"option"]))), sum)

    # gample$family$getTheta() for Plaice in Jurgens GAMs is 1.195788
    catchSigma(sp) <- array(theta, dim=dim(catchMean(sp)),dimnames=dimnames(catchMean(sp)))
    
    return(sp)
    
}

#---------------------------------------------------------------------------------------------------------
# Extract results to a data frame
#---------------------------------------------------------------------------------------------------------

extract_dsvm_res <- function(z, control, ages, season){
  #detach("package:FLCore", unload=TRUE)
  simNumber <-control@simNumber
  sp        <- c("sp1","sp2","sp3","sp4","sp5")
  dsvm_res             <- as.data.frame(rbind(as.matrix(spp1Landings(sim(z))),
                                              as.matrix(spp2Landings(sim(z))),
                                              as.matrix(spp3Landings(sim(z))),
                                              as.matrix(spp4Landings(sim(z))),
                                              as.matrix(spp5Landings(sim(z)))))
  names(dsvm_res)      <- "landings.wt"
  dsvm_res$discards.wt <- c(rbind(as.matrix(spp1Discards(sim(z))),
                                  as.matrix(spp2Discards(sim(z))),
                                  as.matrix(spp3Discards(sim(z))),
                                  as.matrix(spp4Discards(sim(z))),
                                  as.matrix(spp5Discards(sim(z)))))
  dsvm_res$catch.wt    <- dsvm_res$ landings + dsvm_res$discards
  dsvm_res$effort      <- rep(rep(as.matrix(effort(sim(z))),each=length(ages)),length(sp))
  dsvm_res$option      <- rep(rep(as.matrix(choice(sim(z))),each=length(ages)),length(sp))
  dsvm_res$spp         <- as.factor(c(rep(sp, each=(simNumber*length(ages)*length(season)))))
  dsvm_res$cat         <- ages
  dsvm_res$season      <- c(rep(season, each=simNumber*length(ages)))
  dsvm_res$vessel      <- rep(1:simNumber,each=length(ages))
  dsvm_res$option[is.na(dsvm_res$option)] <- "Stay in port"
  dsvm_res[c(1:4)]     <- lapply(dsvm_res[c(1:4)], function(x) as.numeric(as.character(x)))
  dsvm_res[c(5:9)]     <- lapply(dsvm_res[c(5:9)], function(x) as.factor(x))
  is.num               <- sapply(dsvm_res, is.numeric)
  dsvm_res[is.num]     <- lapply(dsvm_res[is.num], round, 6)
  # # Just focus on sp1 and sp2
  # dsvm_res             <- subset(dsvm_res,(spp %in% c("sp1", "sp2")))
  # trip                 <- count(dsvm_res,c("spp","cat","season","option"))
  # names(trip)[5]       <- "trip"
  # dsvm_res             <- aggregate(cbind(landings.wt, discards.wt, catch.wt, effort)~ spp+cat+season+option, FUN=sum, data=dsvm_res)
  # dsvm_res             <- merge(dsvm_res, trip, by=c("spp","cat", "season","option"),all.x=TRUE)
  return(dsvm_res)
}
