library(RDynState5NAI)
hke <-  mac <-  meg <-  hom <- mon <- new("DynStateInput")

catchMean(hke)  <- array(c(rep(10,2), rep(20,2),rep(4,4)), dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))
catchSigma(hke) <- array(0.2                    ,dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))
catchMean(mac)  <- array(0                      ,dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))
catchSigma(mac) <- array(0.1                    ,dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))
catchMean(meg)  <- array(0                      ,dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))
catchSigma(meg) <- array(0.1                    ,dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))
catchMean(hom)  <- array(0                      ,dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))
catchSigma(hom) <- array(0.1                    ,dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))
catchMean(mon)  <- array(0                      ,dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))
catchSigma(mon) <- array(0.1                    ,dim=c(2,2,2), dimnames=list(cat=1:2,season=as.character(1:2),option=c("a","b")))


effort <- array(c(rep(8,4)), dim=c(2,2), dimnames=list(option=c("a","b"), seasons=as.character(1:2)))

hkePrice <- array(c(8,25),dim=c(2,2), dimnames=list(cat=1:2, seasons=as.character(1:2)))
macPrice <- array(c(0,0),dim=c(2,2), dimnames=list(cat=1:2, seasons=as.character(1:2)))
megPrice <- array(c(0,0),dim=c(2,2), dimnames=list(cat=1:2, seasons=as.character(1:2)))
homPrice <- array(c(0,0),dim=c(2,2), dimnames=list(cat=1:2, seasons=as.character(1:2)))
monPrice <- array(c(0,0),dim=c(2,2), dimnames=list(cat=1:2, seasons=as.character(1:2)))


##############################################################################
# test 1 unconstrained hake quotum, high hake price
###############################################################################

control <- DynState.control(spp1LndQuota= 60,  spp2LndQuota= 5,#tn
                            spp1DisQuota= 100, spp2DisQuota= 5, # tn
                            spp1LndQuotaFine= 1E9, spp2LndQuotaFine= 1E9, # tn
                            spp1DisQuotaFine= 0, spp2DisQuotaFine= 0, # tn
                            fuelUse = 1, 
                            fuelPrice = 12, # euros/day
                            gearMaintenance= 2, # 287 euros/day 
                            landingCosts= 1, # euros/tn
                            addNoFishing= TRUE, increments= 6, 
                            interpolationdistance=2,
		            simNumber= 20, numThreads= 24)


z <- DynState(hke, mac, meg, hom, mon, hkePrice, macPrice, megPrice, homPrice, monPrice, effort, control) 


save(z, file="~/Dropbox/BoB/DynState/July_2016/RDynState5NA_Jul2016/script/prueba_Jul.RData")


##############################################################################
# Check results
###############################################################################
layout (matrix(c(1, 4, 1, 4, 1, 5, 1, 5, 2, 6, 2, 6, 2, 7, 3, 7, 3, 8, 3, 8), 10, 2, byrow = TRUE))
plotChoice(hke,z)
hist(netRev(z),20,  main="Net Revenue",xlab = "", xlim= c(0,900))
plot(apply(effort(sim(z))[,,],1,cumsum)[,1], type="s", ylim=c(0,50), main="Effort",xlab = "season",ylab = "")
for (ii in 2:dim(effort(sim(z)))[2]) lines(apply(effort(sim(z))[,,],1,cumsum)[,ii], type="s", ylim=c(0,50))
hist(apply(spp1Landings(sim(z)),2,sum),20, col="grey", main="Landings HKE",xlab = "", xlim= c(0,60))
hist(apply(spp2Landings(sim(z)),2,sum),20, col="grey", main="Landings MAC",xlab = "", xlim= c(0,60))
hist(apply(spp3Landings(sim(z)),2,sum),20, col="grey", main="Landings MEG",xlab = "", xlim= c(0,60))
hist(apply(spp4Landings(sim(z)),2,sum),20, col="grey", main="Landings HOM",xlab = "", xlim= c(0,60))
hist(apply(spp5Landings(sim(z)),2,sum),20, col="grey", main="Landings MON",xlab = "", xlim= c(0,60))

z@Sim@Derivative[,,]
