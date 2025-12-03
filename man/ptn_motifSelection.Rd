\name{ptn_motifSelection}
\alias{ptn_motifSelection,postNetData-method}
\title{Accessor for the "motifSelection" slot of a \code{postNetData} object
}

\description{Retrieve the "motifSelection" slot from a \code{postNetData} object, which holds de-novo identified motifs in the gene sets of interest that have enrichment p-values below the given threshold. 
}

\usage{
\S4method{ptn_motifSelection}{postNetData}(x, region)
}

\arguments{
  \item{x}{a \code{postNetData} object.
  }
  \item{region}{A string specifying the sequence region of interest. This can be one of \code{"UTR5"},
  \code{"UTR3"}, or \code{"CDS"}.
  }
}

\value{The \code{ptn_motifSelection} function returns a character vector of de-novo identified motifs in the gene sets of interest that have enrichment p-values below the threshold set with the \code{stremeThreshold} parameter of the \link[postNet]{motifAnalysis} function. 
}

\references{
Timothy L. Bailey, James Johnson, Charles E. Grant, William S. Noble, "The MEME Suite", Nucleic Acids Research, 43(W1):W39-W49, 2015.

Timothy L. Bailey, "STREME: Accurate and versatile sequence motif discovery", Bioinformatics, 2021. https://doi.org/10.1093/bioinformatics/btab203
}

\seealso{
\link[memes]{memes} \cr
\link[memes]{runStreme} \cr
\link[postNet]{motifAnalysis}
}

\examples{
# Note that as users must provide the path to the MEME Suite executables, it is not possile to
# create a fully executable example as this path will vary depending on many factors. An example of
# how to run this function is provided below, and can be updated with the correct memePath argument. 

\dontrun{
# load example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn
 
 ptn <- motifAnalysis(ptn = ptn, 
                     stremeThreshold = 0.05,
                     minwidth = 6,
                     memePath = "/meme/bin", 
                     region = c('UTR5'),
                     subregion = NULL,
                     subregionSel = NULL)
                     
# Extract the significantly enriched motifs:
 myMotifs <- ptn_motifSelection(ptn = ptn,
                                region = 'UTR5')
                                
 str(myMotifs)
  }
}
