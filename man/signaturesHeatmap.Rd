\name{signaturesHeatmap}
\alias{signaturesHeatmap}
\title{Generate heatmaps to evaluate the regulation of gene signatures in a \code{postNetData} object
}

\description{
}

\usage{
signaturesHeatmap(ptn, 
                  signatureList, 
                  unit = "FDR", 
                  addEff = NULL, 
                  pdfName = NULL)
}

\arguments{
  \item{ptn}{A \code{postNetData} object.
  }
  \item{signatureList}{A named list of vectors of gene IDs   
  }
  \item{unit}{"'unit' can be only 'FDR' or p25, p50, p75 or any percentile in format p with number). 
  If FDR: -log10 FDR wilcoxon test corrected for multitesting with direction of ecdf shift, 
  if percentile: difference from background for one of the percentiles from ecdf"
  }
  \item{addEff}{addEff must be a named list of named numeric vectors.
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
