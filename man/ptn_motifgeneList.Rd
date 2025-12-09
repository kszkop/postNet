\docType{methods}
\name{ptn_motifGeneList}
\alias{ptn_motifGeneList,postNetData-method}
\title{Accessor for the de novo motif analysis results for each gene set of interest in a \code{postNetData} object
}

\description{Retrieve the results of de novo motif identification for each gene set of interest in a \code{postNetData} object. Data for all motifs identified for each gene set of interest in the input \code{geneList} with enrichment p-values below the given threshold are stored in \link[universalmotif]{universalmotif} format. 
}

\usage{
\S4method{ptn_motifGeneList}{postNetData}(x, region, geneList)
}
\arguments{
  \item{x}{A \linkS4class{postNetData} object.
  }
  \item{region}{A string specifying the sequence region of interest. This can be one of \code{"UTR5"},
  \code{"UTR3"}, or \code{"CDS"}.
  }
  \item{geneList}{A string providing the name of the gene set of interest from \code{geneList} to extract de novo
  motif results for. 
  }
}

\value{The \code{ptn_motifGeneList} function returns additional information on the selected de novo identified motifs for the specified gene list in \link[universalmotif]{universalmotif} format.
}

\references{
Timothy L. Bailey, James Johnson, Charles E. Grant, William S. Noble, "The MEME Suite", Nucleic Acids Research, 43(W1):W39-W49, 2015.

Timothy L. Bailey, "STREME: Accurate and versatile sequence motif discovery", Bioinformatics, 2021. https://doi.org/10.1093/bioinformatics/btab203
}

\seealso{
\link[memes]{memes} \cr
\link[memes]{runStreme} \cr
\link[universalmotif]{universalmotif} \cr
\link[postNet]{motifAnalysis}
}

\examples{
# Note that as users must provide the path to the MEME Suite executables, it is not possible to
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
                     

# Extract full motif results from one gene set of interest:                                
  motifsTransUp <- ptn_motifGeneList(ptn = ptn,
                                      region = 'UTR5',
                                      geneList = 'translationUp')
                                
 str(motifsTransUp)
  }
}

