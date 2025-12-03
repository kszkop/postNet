\docType{methods}
\name{ptn_codonSelection}
\alias{ptn_codonSelection,postNetData-method}
\title{Accessor for the "codonSelection" slot of a \code{postNetData} object
}

\description{
}

\usage{
ptn_codonSelection(x, comparison)
}

\arguments{
  \item{x}{A \code{postNetData} object.
}
  \item{comparison}{An integer specifying which comparison to select results from. For example, \code{1} would
  select the first comparison specified in the \code{comparisons} parameter of the \link[postNet]{codonUsage} function.
  }
}

\value{The value returned is a named list of the selected codons or AAs that were signficantly enriched or depleted (according to the thresholds selected in the \link[postNet]{codonUsage} function) in the specified comparisons. This list can be used as input for the \link[postNet]{codonCalc} function. 
}

\seealso{
\link[postNet]{codonUsage} \cr
\link[postNet]{codonCalc}
}

\examples{
# load example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn
 
# Run codon usage analysis of single codons
 ptn <- codonUsage(ptn = ptn, 
          annotType = 'ptnCDS', 
          sourceSeq = "load", 
          analysis = 'codon', 
          codonN = 1, 
          pAdj = 0.01, 
          rem5 = TRUE, 
          plotHeatmap = FALSE, 
          thresOddsUp = 0.4, 
          thresFreqUp = 0.4, 
          thresOddsDown = 0.4, 
          thresFreqDown = 0.4, 
          subregion = NULL, 
          subregionSel = NULL, 
          comparisons = list(c(1,2)), 
          plotType_index = "violin", 
          pdfName = 'codons')
          
# Select codons of interest that were significantly enriched or depleted with high frequency, and the highest and lowest Odds ratios
 
 codons <- ptn_codonSelection(ptn, comparison = 1)
}
