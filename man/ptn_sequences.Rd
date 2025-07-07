\docType{methods}
\name{ptn_sequences}
\alias{ptn_sequences,postNetData-method}
\title{Accessor for the reference annotation sequences in the sequence annotation slots of a \code{postNetData} object
}

\description{The "postNetAnnot" slot holds the reference sequence annotations as \code{postNetRegion} lists. Each annotation for the different sequence regions (5'UTR, CDS, 3'UTR, and optionally CCDS) are lists of three character vectors ("id", "geneID", and "sequence"). The \code{ptn_sequences} function retrieves sequences from the specified mRNA region. 
}

\usage{
\S4method{ptn_sequences}{postNetData}(x, region)
}

\arguments{
  \item{x}{a \code{postNetData} object.
  }

  \item{region}{The sequence region to be accessed. Can be either: \code{"UTR5"}, \code{"CDS"}, or \code{"UTR3"}.
  }
}

\value{
A character vector containing the reference annotation sequences corresponding to the specified region.
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

# Access the 5'UTR sequences from the postNetData object
  UTR5seqs <- ptn_sequences(ptn, "UTR5")
  str(UTR5seqs)
  
}

\seealso{
\code{\link[postNet]{postNetStart}} for more information on obtaining reference sequence annotations from different sources. \cr
\code{\link[postNet]{postNetData-class}}
}
