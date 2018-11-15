\name{DynState.control-class}
\docType{class}
\alias{DynState.control-class}

\title{Class DynState.control}
\description{ A class to represent the controls used in the Dynamic State Variable Model}
\section{Objects from the Class}{
Objects can be created by calls to the creator of this class in the form \code{FDynState.control()}.
}
\section{Slots}{
  \describe{
    \item{\code{spp1LndQuota}:}{Object of class \code{"numeric"}, size of the spp1 landing quota }
    \item{\code{spp2LndQuota}:}{Object of class \code{"numeric"}, size of the spp2 landing quota }
    \item{\code{spp1DisQuota}:}{Object of class \code{"numeric"}, size of the spp1 discard quota }
    \item{\code{spp2DisQuota}:}{Object of class \code{"numeric"}, size of the spp2 discard quota }
    \item{\code{spp1LndQuotaFine}:}{Object of class \code{"numeric"}, fine per unit landings for exceeding quota spp1 }
    \item{\code{spp2LndQuotaFine}:}{Object of class \code{"numeric"}, fine per unit landings for exceeding quota spp2 }
    \item{\code{spp1DisQuotaFine}:}{Object of class \code{"numeric"}, fine per unit discards for exceeding quota spp1 }
    \item{\code{spp2DisQuotaFine}:}{Object of class \code{"numeric"}, fine per unit discards for exceeding quota spp2 }
    \item{\code{simNumber}:}{Object of class \code{"integer"}, indicating the number of simulations to do in the forward iterations }
    \item{\code{sigma}:}{Object of class \code{"numeric"}, indicating how far the probabilties concentrate in the optimal location }
    \item{\code{fuelUse}:}{Object of class \code{"numeric"}, fuel use per unit of effort }
    \item{\code{fuelPrice}:}{Object of class \code{"numeric"}, fuel price per unit of fuel }   
    \item{\code{landingCosts}:}{Object of class \code{"numeric"}, costs of one unit of landings }
    \item{\code{gearMaintenance}:}{Object of class \code{"numeric"}, costs of gear maintenance }
    \item{\code{increments}:}{Object of class \code{"integer"}, the number of bins for discretizing the catch in the discrete prob distribution}
    \item{spp1Incs}{the number of bins for discretizing the catch in the discrete prob distribution for spp1}
    \item{spp2Incs}{the number of bins for discretizing the catch in the discrete prob distribution for spp2}
    \item{spp3Incs}{the number of bins for discretizing the catch in the discrete prob distribution for spp3}
    \item{spp4Incs}{the number of bins for discretizing the catch in the discrete prob distribution for spp4}
    \item{spp5Incs}{the number of bins for discretizing the catch in the discrete prob distribution for spp5}
    \item{\code{addNoFishing}:}{Object of class \code{"integer"}, describing whether to add a single patch that reperesents no fishing (no effort, no catch) } 
    \item{\code{spp1DiscardSteps}:}{Object of class \code{"integer"}, number of discardsteps for spp1, this results in a multiplication of the number of patches }
    \item{\code{spp2DiscardSteps}:}{Object of class \code{"integer"}, number of discardsteps for spp2, this results in a multiplication of the number of patches }
    \item{\code{spp3DiscardSteps}:}{Object of class \code{"integer"}, number of discardsteps for spp3, this results in a multiplication of the number of patches }
    \item{\code{spp4DiscardSteps}:}{Object of class \code{"integer"}, number of discardsteps for spp4, this results in a multiplication of the number of patches }
    \item{\code{spp5DiscardSteps}:}{Object of class \code{"integer"}, number of discardsteps for spp5, this results in a multiplication of the number of patches }
    \item{\code{choiceDist}:}{Object of class \code{"integer"}, describing whether to output the probability of landings for each of the choices }
    \item{\code{numThreads}:}{Object of class \code{"integer"}, indicating the number of threads used for backward solving problem }

}
}
\section{Methods}{
}
\author{Jan Jaap Poos}
\examples{
}
\keyword{classes} 
