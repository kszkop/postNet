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
  \item{ads}{An S4 object of class \code{Anota2seqDataSet} (resulting from the \code{\link[anota2seq]{anota2seqRun}} function in the \pkg{anota2seq} package). 
  }
  \item{contrast}{An integer selecting the \pkg{anota2seq} contrast to retrieve regulatory effect measurements from. 
  For more details on contrasts, please see the \pkg{anota2seq} vignette. 
  }
  \item{dataName}{A string specifying the name that will be given to data in the \code{anota2seq} object. 
  }
  \item{effects_names}{
  }
  \item{signatureList}{A named list of vectors containing gene IDs for the signatures of interest to be examined. Note that several
  signatures of translational regulation are provided with the package for both human and mouse. These can be retrieved by XXXXX
  }
  \item{generalName}{
  }
  \item{signature_colours}{
  }
  \item{xlim}{
  }
  \item{scatterXY}{
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
}

\examples{
}
