\name{codonCalc}
\alias{codonCalc}
\title{Calculate and visualize codon or amino acid counts and/or frequency
}

\description{The \code{codonCalc} function calculates the number or frequency of selected codons or amino acids for each gene.
In addition, codon or amino acid content for gene sets of interest can be compared against background, or other genes sets. 
Several options are available for plotting the comparisons.
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
  \item{ptn}{A \code{postNetData} object.
  }
  \item{featsel}{A named list of vectors specifying one or more codons or amino acids to qunatify and/or visualize. 
  For example, for use in the downstream \link[postNet]{featureIntegration}, it is usually desireable to select 
  those with high frequency, and the highest and lowest Odds Ratios identified using the \link[postNet]{codonUsage} 
  function. This list can be easily generated using the S4 method \link[postNet]{ptn_codonSelection}. See example below. 
  }
  \item{analysis}{A string specifying the type of the analysis to perform. The options are \code{"codon"}, or \code{"AA"} 
  to assess amino acid usage. The default is \code{"codon"}. 
  }
  \item{unit}{A string specifying the unit to quantify. The options are \code{"count"} or \code{"freq"} for frequency. 
  The default is \code{"count"}.
  }
  \item{comparisons}{A list of numeric vectors indicating which gene sets in \code{regulation} or 
  \code{geneList} should be used in pair-wise comparisons. The numeric vectors should correspond to the
  indexes of the gene sets, where 0 denotes the background. For example, \code{list(c(0,2),c(1,2))} 
  indicates that the second gene set in \code{regulation} or \code{geneList} will be compared to the
  background (0), and the first and second gene sets will be compared to each other. The default is
  \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the 
  values for the selected \code{analysis} and \code{unit} will be returned, and statistical 
  comparisons between gene sets will not be performed. The default is \code{TRUE}.
  }
  \item{plotType}{If \code{plotOut} is \code{TRUE}, indicate the type of output plot to generate. Options are,
  \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"ecdf"}. 
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{The number or frequency of selected codons or amino acids is calculated for each gene, and a two-sided Wilcoxon Rank Sum test 
is performed to identify significant differences between gene sets of interest, or against background. If \code{plotOut} is \code{TRUE},
all indicated comparisons between different sets of genes will be performed and plotted.  
}

\value{A named list of vectors for selected codons or amino acids with the calculated number of counts or frequency for each gene. The elements in each vector are named with the gene ID. The \code{codonCalc} function can also return PDF files of output plots with statistical comparisons between gene sets of interest.
}

\references{
}

\seealso{
\link[postNet]{codonUsage} \cr
\link[postNet]{ptn_codonSelection} \cr
\link[postNet]{ptn_codonAnalysis} 
}

\examples{
# load example data:
 data(postNetExample)
    
# Initialize the postNetData object with custom gene lists:
 myGenes <- postNetExample$geneList
 myBg <- postNetExample$background
 myEffect <- postNetExample$effect

 ptn <- postNetStart(
   geneList = myGenes,
   geneListcolours = c("#7FB7BE","#DB7F67"),
   customBg = myBg,
   effectMeasure = myEffect,
   selection = "random",
   setSeed = 123,
   source = "load", 
   species = "human"
 )
 
# Calculate and plot ecdfs of codons counts, and compare between gene sets of interest

 codons <- list(Up = c('GAA','GAT), Down = c('GAG','CAG'))

 codonCounts  <- codonCalc(ptn = ptn,
                          analysis='codon', 
                          featsel= codons,
                          unit='count', 
                          comparisons = list(c(1,2)),
                          plotType = 'ecdf')
                      
 str(GAAcounts)
}
