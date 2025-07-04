\name{checkAvailableVersions}
\alias{checkAvailableVersions}
\title{Check the available RefSeq release versions for pre-prepared reference sequence annotations}

\usage{
checkAvailableVersions(species)
}

\arguments{
\item{species}{A string specifying the species of interest. Currently, \code{"human"} and \code{"mouse"} are supported.}
}
\description{
The \code{checkAvailableVersions} function retrieves available release versions of pre-prepared reference sequence annotations based on the NCBI RefSeq database.
}

\value{A character vector of available RefSeq release versions for the specified species. Selected versions can be specified using the \code{version} parameter of the \link[postNet]{postNetStart} function. 
}

\references{
See \url{https://www.ncbi.nlm.nih.gov/refseq/} for RefSeq documentation.
}

\seealso{
\link[postNet]{postNetStart} for details on creating reference sequence annotations.
}

\examples{
# List RefSeq release versions of available annotations in the package
checkAvailableVersions(species = 'human')
}
