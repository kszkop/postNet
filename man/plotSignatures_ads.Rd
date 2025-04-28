\name{plotSignatures_ads}
\alias{plotSignatures_ads}
\title{Analysis of signatures from an anota2seq object
}

\description{
}

\usage{
plotSignatures_ads(ads, 
                  contrast, 
                  dataName, 
                  effects_names = c("Total mRNA Log2FC", "Polysome-associated mRNA Log2FC", "Buffering Log2FC", "Translation Log2FC"), 
                  signatureList, 
                  generalName, 
                  signature_colours, 
                  xlim = c(-2, 2), 
                  scatterXY = NULL, 
                  tableCex = 1, 
                  pdfName = NULL)
}

\arguments{
  \item{ads}{An S4 object of class \code{Anota2seqDataSet} (resulting from the \link[anota2seq]{anota2seqRun} function in the \pkg{anota2seq} package). 
  }
  \item{contrast}{An integer selecting the \pkg{anota2seq} contrast to retrieve regulatory effect measurements from. 
  For more details on contrasts, please see the \pkg{anota2seq} vignette. 
  }
  \item{dataName}{A string specifying the name that will be given to selected contrast from the \code{Anota2seqDataSet} object. For example: \code{"Test vs. DMSO"}. 
  }
  \item{effects_names}{A character vector specifying axis names for the fold change scatterplot and ecdf plots. The default is, 
  \code{c("Total mRNA Log2FC", "Polysome-associated mRNA Log2FC", "Buffering Log2FC", "Translation Log2FC")}. 
  }
  \item{signatureList}{A named list of vectors containing gene IDs for the signatures of interest to be examined. Note that several
  signatures of translational regulation are provided with the package for both human and mouse. These can be retrieved using the
  \link[postNet]{retrieveSignatures} function.
  }
  \item{generalName}{A character vector specifying the names of the signatures to be plotted. 
  }
  \item{signature_colours}{A character vector specifying colours for plotting each signature. Note that \code{signature_colours} must be the same length as
  \code{signatureList}. 
  }
  \item{xlim}{A numeric vector specifying the x-axis limits for the ecdf plots. Default is \code{c(-2,2)}.
  }
  \item{scatterXY}{A numeric value specifying the scale of the x and y axes for the fold change scatterplot. For example, \code{scatterXY = 3}, would indicate that both
  the x and y axis will have a Log2FC range from -3 to 3. The default is \code{NULL}. 
  }
  \item{tableCex}{A numeric value specifying the size of the text in the table of statistical results that appears at the top of the ecdf plots.
  Default is \code{1}. 
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{
}

\value{No value is returned. Graphical outputs are generated as PDF files. 
}

\seealso{
\link[postNet]{signaturesHeatmap} \cr
\link[postNet]{plotSignatures} \cr
link[postNet]{retrieveSignatures}
}

\examples{
}
