\docType{methods}
\name{ptn_effect}
\alias{ptn_effect,postNetData-method}
\title{Accessor for the regulatory effect measurement from the "effect" slot of a \code{postNetData} object
}

\description{Retrieve the regulatory effect measurement from a \code{postNetData} object. If using \code{geneList}, the "effect" slot holds a user-provided named numeric vector corresponding to a custom regulation effect measurement. If using \code{ads}, the values in the "effect" slot will correspond to the selected \code{regulationGen} from \code{ads}. 
}

\usage{
\S4method{ptn_effect}{postNetData}(ptn)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }
}

\value{A numeric vector with names corresponding to the gene/transcript IDs in \code{geneList} or \code{ads}.
}

\examples{

# load and create example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn
 
# Access the regulatory effect measurement stored in the postNetData object:
  myEffect_out <- ptn_effect(ptn)
  str(myEffect_out)
}

\seealso{
\code{\link[postNet]{postNetStart}} for details on selecting or providing the regulatory effect measurement. \cr
\code{\link[postNet]{postNetData-class}}
}
