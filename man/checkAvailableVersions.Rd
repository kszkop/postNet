\name{checkAvailableVersions}
\alias{checkAvailableVersions}
\title{List the available RefSeq release versions for pre-prepared reference sequence annotations}

\description{
The \code{checkAvailableVersions} function lists the available NCBI RefSeq release versions for which reference sequence annotations are provided in the package.
}

\usage{
checkAvailableVersions(species)
}

\arguments{
\item{species}{A string specifying the species of interest. Currently, \code{"human"} and \code{"mouse"} are supported.}
}

\value{The available RefSeq release versions for the specified species will be printed to the console. The desired annotation version to be used in analysis can be specified using the \code{version} parameter of the \link[postNet]{postNetStart} function. 
}

\references{
See the RefSeq documentation at \url{https://www.ncbi.nlm.nih.gov/refseq/}.
}

\seealso{
\link[postNet]{postNetStart} for details on creating reference sequence annotations.
}

\examples{
# List RefSeq release versions of available annotations in the package
 checkAvailableVersions(species = "human")
 checkAvailableVersions(species = "mouse")
}
