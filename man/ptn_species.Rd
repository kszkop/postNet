\docType{methods}
\name{ptn_species}
\alias{ptn_species,postNetData-method}
\title{Accessor for the "species" slot of a \code{postNetData} object
}

\description{Retrieve the "species" slot from a \code{postNetData} object, which holds the source species of the RefSeq sequence annotations. Currently, \code{"human"} and \code{"mouse"} are available if \code{source} is \code{"load"}, \code{"create"}, \code{"createFromSourceFiles"}, or \code{"createFromFasta"}. 
}

\usage{
\S4method{ptn_species}{postNetData}(x)
}

\arguments{
  \item{x}{A \linkS4class{postNetData} object.
  }
}

\value{A character specifying the species of origin for the RefSeq annotations used to generate the reference sequences.
}

\examples{

# load and create example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Access species of origin for the RefSeq annotation used to create the reference sequences
  ptn_species(ptn)

}

\seealso{
\code{\link[postNet]{postNetStart}} for details on available species, and providing custom annotations. \cr
\code{\link[postNet]{postNetData-class}}
}
