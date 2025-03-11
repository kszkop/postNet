\name{contentMotifs}
\alias{contentMotifs}
\title{
}

\description{The number of motifs in the selected region for each gene is calculated. There is also an option to predict the presence of g-quadruplexes using the pqsfinder package https://bioconductor.org/packages/release/bioc/html/pqsfinder.html . If plotOut is TRUE, all indicated comparisons between different sets of genes are performed and plotted. 
}

\usage{
contentMotifs(ptn, 
              motifsIn, 
              seqType = "dna", 
              dist = 1, 
              min_score = 47, 
              unitOut = "number", 
              resid = FALSE, 
              region, 
              subregion = NULL, 
              subregionSel = NULL, 
              comparisons = NULL, 
              pdfName = NULL, 
              plotOut = TRUE)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
  \item{motifsIn}{motif sequences to be analyzed. IUPAC code or [ ] can be used for ambiguities. Provide "G4" for 
  g-quadruplexes. G-quadruplexes presence is predicted using pqsfinder package 
  https://bioconductor.org/packages/release/bioc/html/pqsfinder.html 
  }
  \item{seqType}{A string speicfying the type of sequence being provided to search for motifs. The options are,
  \code{"dna"}, \code{"rna"}, or \code{"protein"}. The default is \code{"dna"}.
  }
  \item{dist}{Minimal distance between motifs. The default is 1.
  }
  \item{min_score}{Threshold for the quality of the g-quadruplexes prediction. It is required for "G4" option. 
  The default is 47.
  }
  \item{unitOut}{A string to specify whether the output should be the \code{"number"} of uORFs detected per gene, or
  the \code{"position"} of the detected uORFs in the mRNA sequence (nucleotide position of the start and stop). Note
  that only the output with \code{"number"} can be used in the downstream analysis with the
  \link[postNet]{featureIntegration} function. The default is \code{"number"}.
  }
  \item{resid}{Residuals after correction for the length of the sequences should be used instead. It can be TRUE 
  only if unitOut is "number". The default is FALSE.
  }
  \item{region}{A character vector specifying the sequence region(s) of interest. This can be \code{"UTR5"},
  \code{"UTR3"}, \code{"CDS"}, or any combination of the regions.
  }
  \item{subregion}{Optionally, it is possible to specify a more specific subregion of the sequences to either
  select or exclude from the analysis. This is done by providing a numeric value indicating the number of
  nucleotides from either the start of the sequence region if positive, or the end if negative. For example,
  \code{subregion = 50}, would denote the first 50 nucleotides of the sequence region(s) specified by \code{region},
  while \code{subregion = -50}, would denote the last 50 nucleotides of the sequence region(s). The
  default is \code{NULL}.
  }
  \item{subregionSel}{If a value is provided for \code{subregion}, indicate whether the subregion should be
  selected (limiting the analysis to this region) or excluded by specifying either \code{"select"} or 
  \code{"exclude"}. The default is \code{NULL}.
  }
  \item{comparisons}{A list of numeric vectors indicating which gene sets in \code{regulation} or 
  \code{geneList} should be used in pair-wise comparisons. The numeric vectors should correspond to the
  indexes of the gene sets, where 0 denotes the background. For example, \code{list(c(0,2),c(1,2))} 
  indicates that the second gene set in \code{regulation} or \code{geneList} will be compared to the
  background (0), and the first and second gene sets will be compared to each other. The default is
  \code{NULL}.
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the 
  values specified with \code{unitOut} will be returned, and statistical comparisons between gene sets will not be
  performed. Note that plots will only be generated if \code{unitOut = "position"}. The default is \code{TRUE}. 
  }
}

\details{
}

\value{List of vectors for selected regions with calculated number of motifs present for each gene. The function will also output plots with desired comparisons.
}

\references{
}

\author{
}

\note{
}

\seealso{
}

\examples{

}
