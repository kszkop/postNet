\docType{methods}
\name{ptn_selection}
\alias{ptn_selection,postNetData-method}
\title{Accessor for the isoform selection parameter from the "selection" slot of a \code{postNetData} object
}

\description{Retrieve the "selection" slot from a \code{postNetData} object, which holds the method chosen to select which mRNA isoforms are used in analyses. These can be one of: \code{"longest"}, \code{"shortest"}, or \code{"random"}.
}

\usage{
\S4method{ptn_selection}{postNetData}(x)
}

\arguments{
  \item{x}{A \linkS4class{postNetData} object.
  }
}

\value{A character specifying the method used to select mRNA isoforms. Can be either \code{"longest"}, \code{"shortest"}, or \code{"random"}.
}

\examples{

# load and create example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Access mRNA isoform selection method:
  ptn_selection(ptn)

}

\seealso{
\code{\link[postNet]{postNetStart}} for details on the different mRNA isoform selection options. \cr
\code{\link[postNet]{postNetData-class}}
}

