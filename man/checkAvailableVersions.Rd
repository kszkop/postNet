\name{checkAvailableVersions}
\alias{checkAvailableVersions}
\title{Check the available RefSeq release versions for pre-prepared reference sequence annotations}
\usage{
checkAvailableVersions(
  species = NULL,
)
}
\arguments{
\item{species}{Select the species of interest. Currently, \code{'human'} and \code{'mouse'} are supported.}
}
\description{
The \code{checkAvailableVersions} function retrieves and lists the available release versions for the pre-prepared reference sequence annotations based on the NCBI RefSeq database.
}

\examples{

# See RefSeq release versions of available annotations in the package
checkAvailableVersions(species = 'human')
  
}

\seealso{
\code{\link[postNet]{postNetStart}} for details on creating reference sequence annotations.
}
