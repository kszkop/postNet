\docType{methods}
\name{ptn_colours}
\alias{ptn_colours,postNetData-method}
\title{Accessor for the "colours" slot of a \code{postNetData} object
}

\description{Retrieve the "colours" slot of a a \code{postNetData} object, which holds a character vector with the user-provided colours used for plotting. The vector must be the same length as number of elements in \code{geneList} (e.g., one colour for each gene set of interest).  
}

\usage{
\S4method{ptn_colours}{postNetData}(x)
}

\arguments{
  \item{x}{a \code{postNetData} object.
  }
}

\value{A character vector specifying the colours to be used for plotting.
}

\examples{

# load and create example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Access version of the RefSeq annotation used to create the reference sequences
  myCols <- ptn_colours(ptn)
  str(myCols)

}

\seealso{
\code{\link[postNet]{postNetStart}} for details on selecting colours. \cr
\code{\link[postNet]{postNetData-class}}
}
