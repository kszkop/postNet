\name{codonCalc}
\alias{codonCalc}
\title{Quantify and compare codon or amino acid usage across gene sets
}

\description{After running \code{codonUsage}, the \code{codonCalc} function can be used to calculate the number or frequency of selected codons or amino acids (AA) for each gene. These values can be used in the downstream \link[postNet]{featureIntegration} analysis. In addition, codon or amino acid content for gene sets of interest can be compared against background, or other gene sets. Several options are available for plotting the comparisons.
}

\usage{
codonCalc(ptn, 
          featsel, 
          analysis = "codon", 
          unit = "count", 
          comparisons = NULL, 
          plotOut = TRUE, 
          plotType = "ecdf", 
          pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }
  \item{featsel}{A named list of vectors specifying one or more codons or amino acids to quantify and/or visualize. For example, for use in the downstream \link[postNet]{featureIntegration}, it is usually desirable to select those with high frequency, and the highest and lowest odds ratios identified using the \link[postNet]{codonUsage} function. This list can be easily generated using the S4 method \link[postNet]{ptn_codonSelection}. 
  }
  \item{analysis}{A string specifying whether to assess \code{"codon"}, or \code{"AA"} (amino acid) usage. The default is \code{"codon"}. 
  }
  \item{unit}{A string specifying the unit to quantify. The options are \code{"count"} or \code{"freq"} for frequency. Here, frequency corresponds to the number of codons relative to all codons in the gene. The default is \code{"count"}.
  }
  \item{comparisons}{A list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the values for the selected \code{analysis} and \code{unit} will be returned, and statistical comparisons between gene sets will not be performed. The default is \code{TRUE}.
  }
  \item{plotType}{If \code{plotOut} is \code{TRUE}, indicate the type of output plot to generate. Options are \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"ecdf"}. 
  }
  \item{pdfName}{Name to be appended to output PDF files. The default is \code{NULL}.
  }
}

\details{The number or frequency of selected codons or amino acids is calculated for each gene, and a two-sided Wilcoxon Rank Sum test is performed to identify significant differences between gene sets of interest, or against background. If \code{plotOut} is \code{TRUE}, all indicated comparisons between different sets of genes will be performed and plotted.  
}

\value{A named list of numeric vectors for selected codons or amino acids with the calculated counts or frequencies for each gene. The elements in each vector are named with the gene ID. This list can be used as input for the downstream \link[postNet]{featureIntegration} analysis. The \code{codonCalc} function can also return PDF files of output plots with statistical comparisons between gene sets of interest.
}

\seealso{
\link[postNet]{codonUsage} \cr
\link[postNet]{ptn_codonSelection} \cr
\link[postNet]{ptn_codonAnalysis} 
}

\examples{
# load example data:
 postNetExample <- postNetParcel::postNetExample
 ptn <- postNetExample$ptn
 
# Run codon usage analysis
 ptn <- codonUsage(ptn = ptn, 
          annotType = "ptnCDS", 
          sourceSeq = "load", 
          analysis = "codon", 
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
          pdfName = "codons")
          
# Select codons of interest with high frequency, and the highest and lowest odds ratios
 codons <- ptn_codonSelection(ptn, comparison = 1)
 
# Calculate and plot ecdfs of codons counts, and compare between gene sets of interest
 codonCounts  <- codonCalc(ptn = ptn,
                          analysis = "codon", 
                          featsel = codons,
                          unit = "count", 
                          comparisons = list(c(1,2)),
                          plotType = "ecdf")
                      
 str(codonCounts)
}
