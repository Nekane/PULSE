\name{Sim-class}
\docType{class}
\alias{Sim-class}

\title{Class Sim}
\description{ A class to represent the outputs of the forward iterations from a dynamic state variable model for a mixed fishery under quota constraints. This is part of the returned object from a model run }
\section{Slots}{
  \describe{
    \item{\code{choice}:}{Object of class \code{"array"}, indicating which option is chosen. rows are forward simulation number, colums are timesteps }
    \item{\code{effort}:}{Object of class \code{"array"}, indicating how much effort is used in this time step }
    \item{\code{spp1DisHold}:}{Object of class \code{"array"}, total (annual) spp1 discards after this time step  }
    \item{\code{spp1Discards}:}{Object of class \code{"array"}, the amount of spp1 discards in this time step  }
    \item{\code{spp1Landings}:}{Object of class \code{"array"}, spp1 landings in this time step }
    \item{\code{spp1LndHold}:}{Object of class \code{"array"}, total (annual) spp1 landings after this timestep }
    \item{\code{spp1Rand}:}{Object of class \code{"array"}, random number drawn, in the forward simulation on which the spp1 landings are based }
    \item{\code{spp2DisHold}:}{Object of class \code{"array"}, total (annual) spp2 discards after this time step  }
    \item{\code{spp2Discards}:}{Object of class \code{"array"}, the amount of spp2 discards in this time step  }
    \item{\code{spp2Landings}:}{Object of class \code{"array"}, spp2 landings in this time step }
    \item{\code{spp2LndHold}:}{Object of class \code{"array"}, total (annual) spp2 landings after this timestep }
    \item{\code{spp2Rand}:}{Object of class \code{"array"}, random number drawn, in the forward simulation on which the spp2 landings are based }
    \item{\code{spp3DisHold}:}{Object of class \code{"array"}, total (annual) spp3 discards after this time step  }
    \item{\code{spp3Discards}:}{Object of class \code{"array"}, the amount of spp3 discards in this time step  }
    \item{\code{spp3Landings}:}{Object of class \code{"array"}, spp3 landings in this time step }
    \item{\code{spp3LndHold}:}{Object of class \code{"array"}, total (annual) spp3 landings after this timestep }
    \item{\code{spp3Rand}:}{Object of class \code{"array"}, random number drawn, in the forward simulation on which the spp3 landings are based }
    \item{\code{spp4DisHold}:}{Object of class \code{"array"}, total (annual) spp4 discards after this time step  }
    \item{\code{spp4Discards}:}{Object of class \code{"array"}, the amount of spp4 discards in this time step  }
    \item{\code{spp4Landings}:}{Object of class \code{"array"}, spp4 landings in this time step }
    \item{\code{spp4LndHold}:}{Object of class \code{"array"}, total (annual) spp4 landings after this timestep }
    \item{\code{spp4Rand}:}{Object of class \code{"array"}, random number drawn, in the forward simulation on which the spp4 landings are based }
    \item{\code{spp5DisHold}:}{Object of class \code{"array"}, total (annual) spp5 discards after this time step  }
    \item{\code{spp5Discards}:}{Object of class \code{"array"}, the amount of spp5 discards in this time step  }
    \item{\code{spp5Landings}:}{Object of class \code{"array"}, spp5 landings in this time step }
    \item{\code{spp5LndHold}:}{Object of class \code{"array"}, total (annual) spp5 landings after this timestep }
    \item{\code{spp5Rand}:}{Object of class \code{"array"}, random number drawn, in the forward simulation on which the spp5 landings are based }
  }
}
\section{Methods}{
}
\author{Jan Jaap Poos}
\examples{
}
\keyword{classes} 
