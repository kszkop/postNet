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
  \item{signatureList}{A named list of vectors containing gene IDs for the signatures of interest to be examined. Note that several
  signatures of translational regulation are provided with the package for both human and mouse. These can be retrieved using the
  \link[postNet]{retrieveSignatures} function.
  }
  \item{unit}{A string specifying the unit to be plotted in the heatmap. The options are \code{"FDR"},
  or any percentile in format p with a number (for example \code{p25}, \code{p75}, etc.).
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{When \code{unit = 'FDR'}, the values display in the heatmap are obtained from a two-sided Wilcoxon Rank Sum test comparing 
the regulatory effect measurement values for the gene signature of interest against the background. See \link[postNet]{postNetStart} 
for details on background gene sets. P-values are then corrected for multiple testing using the Benjamini & Hochberg method, and the
-log10 FDR value is multiplied by either 1 or -1 corresponding to up- or down-regulation of the gene signature compared to background,
respectively.

When \code{unit = 'p75'} (or any other percentile value), the values displayed in the heatmap are obtained by taking the difference between
the Empirical Cumulative Distribution Function (ecdf) of the regulatory effect measurements for the gene signature and the background,
at the percentile specified.
}

\value{No value is returned. Graphical outputs are generated as PDF files. 
}

\seealso{
\link[postNet]{plotSignatures} \cr
\link[postNet]{plotSignatures_ads} \cr
}

\examples{
# load example data:
 postNetExample <- postNetParcel::postNetExample
 ptn <- postNetExample$ptn
 
# load signature data:
 signatures <- get_signatures("human")
 
# Examine the regulation of various transcripts sensitive to translational regulation by various pathways and factors 
 signaturesHeatmap(ptn,
                  signatureList = signatures,
                  unit = 'FDR',
                  pdfName = "mRNAtranslationSignatures")
}
