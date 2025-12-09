\name{plotSignatures}
\alias{plotSignatures}
\title{Analysis of regulation of gene signatures in a \code{postNetData} object
}

\description{The \code{plotSignatures} function assesses changes in the regulation of known gene signatures in the data set of interest contained in a \code{postNetData} object. Outputs include statistical analyses and visualizations of signature regulation. 
}

\usage{
plotSignatures(ptn, 
              signatureList, 
              signature_colours = NULL,
              dataName, 
              generalName,
              xlim = c(-2, 2), 
              tableCex = 0.7, 
              pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }
  \item{signatureList}{A named list of vectors containing gene IDs for the signatures of interest to be examined. Note that several signatures of translational regulation are provided with the package for both human and mouse. These can be retrieved using the \link[postNet]{get_signatures} function.
  }
  \item{signature_colours}{An optional character vector specifying colours for plotting each signature. Note that \code{signature_colours} must be the same length as \code{signatureList}. The default is \code{NULL}. 
  }
  \item{dataName}{A character string specifying a name for the data set in which the gene signature(s) will be examined. This name will be appended to the name of the PDF file generated. 
  }
   \item{generalName}{A character string specifying a name for the signatures that will be examined. This name will be appended to the name of the PDF file generated. 
  }
  \item{xlim}{A numeric vector specifying the x-axis limits for the eCDF plots. Default is \code{c(-2, 2)}.
  }
  \item{tableCex}{A numeric value specifying the size of the text in the table of statistical results that appears at the top of the eCDF plots. Default is \code{0.7}.
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{Gene signatures of interest are visualized using empirical cumulative distribution functions (eCDFs) of the regulatory effect measurement in a \code{postNetData} object (often log2 fold changes). ECDFs for genes belonging to the signatures of interest are compared to those for all other genes (background). Differences in the regulatory effect measurement distributions are calculated at the quantiles. Significant directional shifts in the eCDFs for gene signature versus background are identified using a Wilcoxon rank-sum test.
}

\value{No value is returned. Graphical outputs are generated as PDF files. 
}

\seealso{
\link[postNet]{signaturesHeatmap} \cr
\link[postNet]{plotSignatures_ads} \cr
\link[postNet]{get_signatures}
}

\examples{
# load example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn
 
# load signature data:
 signatures <- get_signatures("human")
 
# Examine the regulation of mTOR-sensitive transcripts in the regulatory effects stored in a postNetData object:
 plotSignatures(ptn = ptn,
                signatureList = signatures[c("Gandin_etal_2016_mTOR_transUp", "Gandin_etal_2016_mTOR_transDown")],
                signature_colours = c("red", "blue"), 
                dataName = "Example_data",
                generalName = "mTOR_sensitive_translation",
                xlim = c(-3, 3), 
                tableCex = 1,
                pdfName = "mTORsignatures")
}
