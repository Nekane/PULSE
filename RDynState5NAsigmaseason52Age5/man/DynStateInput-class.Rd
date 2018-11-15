\name{DynStateInput-class}
\docType{class}
\alias{DynStateInput-class}

\title{Class DynStateInput}
\description{ A class to represent expected catchrates with mean and CV  }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DynStateInput", ...)}.
}
\section{Slots}{
  \describe{
   
    \item{\code{catchMean}:}{A 3 dimensional object of class \code{"array"}, with the mean catch rate. the first defines three different categories, the second dimension defines the seasons, and the third dimension defines the different options. In the DynState function the mean and sigma are used in the pnorm function.}
    \item{\code{catchSigma}:}{A 3 dimensional object of class \code{"array"}, with the sigma of the  catch rate. the first defines three different categories, the second dimension defines the seasons, and the third dimension defines the different options. In the DynState function the mean and sigma are used in the pnorm function. }
  }
}
\section{Methods}{

}
\author{Jan Jaap Poos}
\examples{
}
\keyword{classes}
