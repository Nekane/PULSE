% $ Id: $
\name{DynState.control}
\alias{DynState.control}
\title{5 species 2 categories quota dynamic state variable model}
\description{
	The \code{DynState.control} function simulates decision and outcomes in a 3 species mixed  
	fisheries under Individual Tradeable Quota constraints or effort constraints. The simulation
        technique used is called dynamic state variable modelling (Clark and Mangel,2000). The input to the model is 
	defined by two DynStateInput objects.
}
\usage{
	DynState.control(spp1LndQuota, spp2LndQuota, spp1DisQuota, spp2DisQuota, spp1LndQuotaFine,spp2LndQuotaFine, spp1DisQuotaFine,spp2DisQuotaFine, simNumber, fuelUse,fuelPrice, landingCosts, increments,addNoFishing, spp1DiscardSteps,spp2DiscardSteps, spp3DiscardSteps,spp4DiscardSteps, spp5DiscardSteps, choiceDist, numThreads)
}

\arguments{
  \item{spp1LndQuota}{The species 1 landing quota}
  \item{spp2LndQuota}{The species 2 landing quota}
  \item{spp1LndQuotaFine}{The fine per unit landings for exceeding species 1 quota}
  \item{spp2LndQuotaFine}{The fine per unit landings for exceeding species 2 quota}
  \item{simNumber}{The number of forward simulations}
  \item{sigma}{The tuning parameter that measures how important is to be near optimal}
  \item{fuelUse}{the fuel use per unit of effort}
  \item{fuelPrice}{the price for one unit of fuel}
  \item{landingCosts}{the cost of one unit of landing}
  \item{gearMaintenance}{the costs of gear maintenance }
  \item{increments}{the number of bins for discretizing the catch in the discrete prob distribution}
  \item{addNoFishing}{integer [0,1], indicating whether to add a single patch that reflects no fishing, catch=0, effort=0 }
  \item{spp1DiscardSteps}{number of discards option for discarding species 1}
  \item{spp2DiscardSteps}{number of discards option for discarding species 2}
  \item{spp3DiscardSteps}{number of discards option for discarding species 3}
  \item{spp4DiscardSteps}{number of discards option for discarding species 4}
  \item{spp5DiscardSteps}{number of discards option for discarding species 5}
  \item{choiceDist}{ouptut the treholds in the prob distribution for the landings?}
  \item{numThreads}{The number of threads used for backward solving DSVM}
}
\details{
	This function returns an object of class dynState with the decision array and the individual
        simulations through this array with information on the effort, catch rates etc.
   }
\value{
    An object of class DynState.
}
\references{}
\author{Jan Jaap Poos}
\note{}
\seealso{}
\examples{}
\keyword{models} 
