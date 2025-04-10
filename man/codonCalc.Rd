\name{codonCalc}
\alias{codonCalc}
\title{Calculate and visualize codon counts and/or frequency
}

\description{
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
  those with high frequency, and the highest and lowest Odds Ratios identified using the \link[postNet{codonUsage} 
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
  values for the selected \code{analysis} and \code{unit} will be retruned, and statistical 
  comparisons between gene sets will not be performed. The default is \code{TRUE}.
  }
  \item{plotType}{If \code{plotOut} is \code{TRUE}, indicate the type of output plot to generate. Options are,
  \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"ecdf"}. 
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{
}

\value{
}

\references{
}

\seealso{
}

\examples{
}
