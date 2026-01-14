\docType{methods}
\name{ptn_geneList}
\alias{ptn_geneList,postNetData-method}
\title{Accessor for the "geneList" slot of a \code{postNetData} object
}

\description{Retrieve the \code{geneList} slot from a \code{postNetData} object, which holds a user-provided named list of one or more character vectors corresponding to gene sets of interest to be examined in the postNet analysis.
}

\usage{
\S4method{ptn_geneList}{postNetData}(ptn)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
}

\value{A named list of one or more character vectors corresponding to gene sets of interest.
}

\examples{

# load and create example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Access the gene sets of interest stored in the postNetData object:
  myGeneList <- ptn_geneList(ptn)
  str(myGeneList)
}

\seealso{
\code{\link[postNet]{postNetStart}} for details on providing gene sets of interest. \cr
\code{\link[postNet]{postNetData-class}}
}
