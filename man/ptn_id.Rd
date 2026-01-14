\docType{methods}
\name{ptn_id}
\alias{ptn_id,postNetData-method}
\title{Accessor for the transcript IDs in the "id" slot of a \code{postNetData} object
}
\description{Retrieve the transcript IDs for the reference annotation sequences from a \code{postNetData} object. The annotations for the different sequence regions (5'UTR, CDS, 3'UTR, and optionally CCDS) are lists of three character vectors ("id", "geneID", and "sequence"). The \code{ptn_id} function retrieves the transcript IDs for sequences from the specified mRNA region. 
}

\usage{
\S4method{ptn_id}{postNetData}(ptn, region)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }

  \item{region}{The sequence region from which to access the transcript IDs. Can be either: \code{"UTR5"}, \code{"CDS"}, or \code{"UTR3"}.
  }
}

\value{
A character vector containing the transcript IDs corresponding to the reference sequence annotations for the specified mRNA region.
}

\examples{

# load and create example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Access the transcript IDs corresponding to the reference 5'UTR sequences from the postNetData object
  myTransIDs <- ptn_id(ptn, "UTR5")
  str(myTransIDs)
  
}

\seealso{
\code{\link[postNet]{postNetStart}} for more information on reference sequence annotations. \cr
\code{\link[postNet]{postNetData-class}}
}
