\docType{methods}
\name{ptn_geneID}
\alias{ptn_geneID,postNetData-method}
\title{Accessor for the gene IDs in the "geneID" slot of a \code{postNetData} object
}
\description{Retrieve the gene IDs for the reference annotation sequences from a \code{postNetData} object. Each annotation for the different sequence regions (5'UTR, CDS, 3'UTR, and optionally CCDS) are lists of three character vectors ("id", "geneID", and "sequence"). The \code{ptn_geneID} function retrieves the gene IDs for sequences from the specified mRNA region. 
}
\usage{
\S4method{ptn_geneID}{postNetData}(x, region)
}

\arguments{
  \item{x}{a \code{postNetData} object.
  }

  \item{region}{The sequence region from which to access the transcript IDs. Can be either: \code{"UTR5"}, \code{"CDS"}, or
  \code{"UTR3"}.
  }
}

\value{
A character vector containing the gene IDs corresponding to the reference sequence annotations for the specified
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
  myGeneIDs <- ptn_geneID(ptn, "CDS")
  str(myGeneIDs)
  
}

\seealso{
\code{\link[postNet]{postNetStart}} for more information on reference sequence annotations. \cr
\code{\link[postNet]{postNetData-class}}
}