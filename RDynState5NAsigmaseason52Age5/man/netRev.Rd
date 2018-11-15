% $ Id: $
\name{netRev}
\alias{netRev}

\title{calculate net revenue from 3 species size structured quota dynamic state variable model}
\description{

	The \code{netRev} function returns net revenue from a dynstate object 

}

\usage{	netRev(dynstate)

}

\arguments{

  \item{dynstate}{A DynState object.}

}

\details{

	This function returns an object of class numeric with the calculated net revenue from a dynstate model. Net revenue is defined as gross revenue - (effort * fueluse * fuelprice) - (landings * landingsCosts) - annual fine.   

   }

\value{

    An object of class numeric.

}

\references{}

\author{Jan Jaap Poos}

\note{}

\seealso{}

\examples{}

\keyword{models} 
