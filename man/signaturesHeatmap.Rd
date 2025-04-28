\name{signaturesHeatmap}
\alias{signaturesHeatmap}
\title{Generate heatmaps to evaluate the regulation of gene signatures in a \code{postNetData} object
}

\description{The \code{signaturesHeatmap} function produces a heatmap of the regulation of gene signatures of interest in a 
\code{postNetData} object. Regulation of a gene signature of interest in the data set can be evaluated using the FDR of a Wilcoxon
Rank Sum test comparing the regulatory effect measurement for the gene signature against the background. Alternatively, the difference
in the regulatory effect measurement for the genes in the signature compared to background can also be visualized at percentiles.  
}

\usage{
signaturesHeatmap(ptn, 
                  signatureList, 
                  unit = "FDR", 
                  pdfName = NULL)
}

\arguments{
  \item{ptn}{A \code{postNetData} object.
  }
  \item{signatureList}{A named list of vectors containing gene IDs for the signatures of interest to be examined.    
  }
  \item{unit}{A string specifying the unit to be plotted in the heatmap. The options are \code{"FDR"},
  or any percentile in format p with a number (for example \code{p25}, \code{p75}, etc.).
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{
}

\value{
}

\seealso{
}

\examples{
}
