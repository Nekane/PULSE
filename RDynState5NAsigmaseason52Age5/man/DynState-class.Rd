\name{DynState-class}
\docType{class}
\alias{DynState-class}

\title{Class DynState}
\description{ A class to represent the different results from the forward and backward calculations of a Dynamic State Variable Model}
\section{Objects from the Class}{
This object is created by the RDynState function
}
\section{Slots}{
  \describe{
    \item{\code{sim}:}{Object of class \code{"sim"}, the result of the forward iterations } 
    \item{\code{control}:}{Object of class \code{"DynState.control"} being the control object that was used to run the model} 
    \item{\code{choiceDist}:}{Object of class \code{"array"}, giving the probability distribution of the landings for the different options }
    \item{\code{spp1Price}:}{Object of class \code{"array"}, the prices that were used as input for species 1}
    \item{\code{spp2Price}:}{Object of class \code{"array"}, the prices that were used as input for species 2 }
    \item{\code{spp3Price}:}{Object of class \code{"array"}, the prices that were used as input for species 3 }
    \item{\code{spp4Price}:}{Object of class \code{"array"}, the prices that were used as input for species 4 }
    \item{\code{spp5Price}:}{Object of class \code{"array"}, the prices that were used as input for species 5 }
  }
}
\section{Methods}{
}
\author{Jan Jaap Poos}
\examples{
}
\keyword{classes} 
