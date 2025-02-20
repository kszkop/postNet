\docType{methods}
\name{ptn_id}
\alias{ptn_id,postNetData-method}
\title{Accessor for the transcript IDs in the "id" slot of a \code{postNetData} object
}
\description{Retrieve the transcript IDs for the reference annotation sequences from a \code{postNetData} object. Each annotation for the different sequence regions (5'UTR, CDS, 3'UTR, and optionally CCDS) are lists of three character vectors ("id", "geneID", and "sequence"). The \code{ptn_id} function retrieves the transcript IDs for sequences from the specified mRNA region. 
}

\usage{
\S4method{ptn_id}{postNetData}(x, region)
}

\arguments{
  \item{x}{a \code{postNetData} object.
  }

  \item{region}{The sequence region from which to access the transcript IDs. Can be either: \code{"UTR5"}, \code{"CDS"}, or
  \code{"UTR3"}.
  }
}

\value{
A character vector containing the transcript IDs corresponding to the reference sequence annotations for the specified
mRNA region.
}

\examples{

# load and create example data:
 data(postNetExample)

 myGenes <- postNetExample$geneList
 myBg <- postNetExample$background
 myEffect <- postNetExample$effect

# Initialize the postNetData object with custom gene lists:
 ptn <- postNetStart(
   geneList = myGenes,
   geneListcolours = c("#7FB7BE","#DB7F67"),
   customBg = myBg,
   effectMeasure = myEffect,
   selection = "random",
   source = "load", 
   species = "human"
 )

# Access the transcript IDs corresonding to the reference 5'UTR sequences from the postNetData object
  myTransIDs <- ptn_id(ptn, "UTR5")
  str(myTransIDs)
  
}

\seealso{
\code{\link[postNet]{postNetStart}} for more information on reference sequence annotations. \cr
\code{\link[postNet]{postNetData-class}}
}
