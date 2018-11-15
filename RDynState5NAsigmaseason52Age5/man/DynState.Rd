% $ Id: $
\name{DynState}
\alias{DynState}

\title{5 species size structured quota dynamic state variable model}
\description{

	The \code{DynState} function simulates decision and outcomes in a 5 species length structured mixed fisheries under Individual Tradeable Quota constraints or effort constraints. The simulation technique used is called dynamic state variable modelling (Clark and Mangel,2000). The input to the model is defined by three DynStateInput objects for statistical distributions of catches for the five species. Prices per category are given in arrays 

}

\usage{	DynState(inputSpp1,inputSpp2, inputSpp3, inputSpp4, inputSpp5, spp1Price, spp2Price, spp3Price, spp4Price, spp5Price, inputEffort, control)

}

\arguments{

  \item{inputSpp1}{A DynStateInput object holding the expected mean catch rate and sd for species 1 (arrays).}

  \item{inputSpp2}{A DynStateInput object holding the expected mean catch rate and sd for species 2 (arrays).}
 
  \item{inputSpp3}{A DynStateInput object holding the expected mean catch rate and sd for species 3 (arrays).}
  
  \item{inputSpp4}{A DynStateInput object holding the expected mean catch rate and sd for species 4 (arrays).}
 
  \item{inputSpp5}{A DynStateInput object holding the expected mean catch rate and sd for species 5 (arrays).}

 \item{spp1Price}{An array holding the prices, structured by category for species 1 (arrays).}

 \item{spp2Price}{An array holding the prices, structured by category for species 2 (arrays).}

 \item{spp3Price}{An array holding the prices, structured by category for species 3 (arrays).}

 \item{spp4Price}{An array holding the prices, structured by category for species 4 (arrays).}

 \item{spp5Price}{An array holding the prices, structured by category for species 5 (arrays).}

  \item{inputEffort}{A DynStateInput object holding the effort needed (array).}

  \item{control}{A control object of class DynState.control}

}

\details{

	This function returns an object of class DynState with the decision array and the individual

        simulations through this array with information on the effort, catch rates and exceeding of quota.

   }

\value{

    An object of class DynState.

}

\references{}

\author{Jan Jaap Poos}

\note{}

\seealso{}

\examples{

hke <-  mac <-  meg <-  hom <- mon <- new("DynStateInput")

catchMean(hke)  <- array(c(rep(12,4),rep(10,4)), dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))
catchSigma(hke) <- array(1,                      dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))
catchMean(mac)  <- array(c(rep(10,4),rep(12,4)), dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))
catchSigma(mac) <- array(1,                      dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))
catchMean(meg)  <- array(0,                      dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))
catchSigma(meg) <- array(0.1,                    dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))
catchMean(hom)  <- array(0,                      dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))
catchSigma(hom) <- array(0.1,                    dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))
catchMean(mon)  <- array(0,                      dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))
catchSigma(mon) <- array(0.1,                    dim=c(2,2,2),dimnames=list(cat=1:2,season=as.character(1:2),option=c("VIIIabd","VII")))

effort <- array(c(rep(8,4)), dim=c(2,2), dimnames=list(option=c("VIIIabd","VII"), seasons=as.character(1:2)))

hkePrice <- array(c(8,8),dim=c(2,2), dimnames=list(cat=1:2,seasons=as.character(1:2)))
macPrice <- array(c(7,7),dim=c(2,2), dimnames=list(cat=1:2,seasons=as.character(1:2)))
megPrice <- array(c(0,0),dim=c(2,2), dimnames=list(cat=1:2,seasons=as.character(1:2)))
homPrice <- array(c(0,0),dim=c(2,2), dimnames=list(cat=1:2,seasons=as.character(1:2)))
monPrice <- array(c(0,0),dim=c(2,2), dimnames=list(cat=1:2,seasons=as.character(1:2)))

control <- DynState.control(spp1LndQuota= 200,  spp2LndQuota= 200,#tn
                            spp1LndQuotaFine= 1E9, spp2LndQuotaFine= 1E9, # tn
                            fuelUse = 1, 
                            fuelPrice = 124, # euros/day
                            landingCosts= 121, # euros/tn
                            gearMaintenance= 271,
                            addNoFishing= TRUE, increments= 6,
                            sigma= 0.1,
                            spp1DiscardSteps=0, spp2DiscardSteps=0,
                            simNumber=20, numThreads=24)


z <- DynState(hke, mac, meg, hom, mon, hkePrice, macPrice, megPrice, homPrice, monPrice, effort, control)



}

\keyword{models} 
