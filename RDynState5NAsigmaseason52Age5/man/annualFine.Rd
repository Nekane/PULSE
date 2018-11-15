% $ Id: $
\name{annualFine}
\alias{annualFine}

\title{calculate total fine from dynamic state variable model }
\description{

	The \code{annualFine} function returns total fine from dynamic state variable model based on overshooting quota for spp1 and fishing effort 

}

\usage{	annualFine(dynstate)

}

\arguments{

  \item{dynstate}{A DynState object.}

}

\details{

	This function returns an object of class numeric with the calculated total fine from a dynstate model. The fine is calculated as the amount of overshooting * Fine.   

   }

\value{

    An object of class numeric.

}

\references{}

\author{Jan Jaap Poos}

\note{}

\seealso{}

\examples{


}

\keyword{models} 
