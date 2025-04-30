\name{plotSignatures_ads}
\alias{plotSignatures_ads}
\title{Analysis of regulation of gene signatures from an anota2seq object
}

\description{The \code{plotSignatures_ads} function assesses changes in the regulation of known gene signatures of interest in an 
\code{Anota2seqDataSet} object from an \pkg{anota2seq} analysis.   
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

\details{Gene signatures of interest are visualized in a scatterplot of log2 fold changes in polysome-associated and total mRNA. The empirical 
cumulative distribution functions (ECDFs) of log2 fold changes in polysome-associated mRNA, total mRNA, translation, and buffering (offsetting) 
are also plotted independently. Fold change ECDFs for genes belonging to the signatures of interest are compare to those for all other genes 
(background). Differences in the fold change distributions are calculated at the quantiles. Significant directional shifts in the ECDFs for gene 
signatures versus background are identified using a Wilcoxon rank-sum test.
}

\value{No value is returned. Graphical outputs are generated as PDF files. 
}

\seealso{
\link[anota2seq]{anota2seqRun} \cr
\link[postNet]{signaturesHeatmap} \cr
\link[postNet]{plotSignatures} \cr
link[postNet]{retrieveSignatures}
}

\examples{
}
