\name{plotSignatures}
\alias{plotSignatures}
\title{Analysis of signatures from a postnetdata object
}

\description{
}

\usage{
plotSignatures(ptn, 
              signatureList, 
              signature_colours = NULL, 
              xlim = c(-2, 2), 
              tableCex = 1, 
              pdfName = NULL)
}

\arguments{
  \item{ptn}{A \code{postNetData} object
  }
  \item{signatureList}{A named list of vectors containing gene IDs for the signatures of interest to be examined. Note that several
  signatures of translational regulation are provided with the package for both human and mouse. These can be retrieved using the
  \link[postNet]{retrieveSignatures} function.
  }
  \item{signature_colours}{A character vector specifying colours for plotting each signature. Note that \code{signature_colours} must be the same length as
  \code{signatureList}. 
  }
  \item{xlim}{A numeric vector specifying the x-axis limits for the ecdf plots. Default is \code{c(-2,2)}.
  }
  \item{tableCex}{
  }
  \item{pdfName}{
  }
}

\details{
}

\value{No value is returned. Graphical outputs are generated as PDF files. 
}

\seealso{
\link[postNet]{signaturesHeatmap} \cr
\link[postNet]{plotSignatures_ads} \cr
link[postNet]{retrieveSignatures}
}

\examples{
}
