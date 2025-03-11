\name{contentMotifs}
\alias{contentMotifs}
\title{Detect and quantify the presence of sequence motifs
}

\description{The \code{contentMotifs} function detects and quantifies the presence of motifs in sequence regions of interest. Known motifs can be provided directly, or the output of the \link[postNet]{motifAnalysis} function, which detects de-novo motifs using STREME (part of the MEME-Suite), can be supplied. Finally, the function can also apply \link[pqsfinder]{pqsfinder} to indeitify G-quadruplexes. In addition, the motifs can be plotted and statistical comparisons can be made between gene sets of interest against background, and/or other genes sets.
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
  \item{motifsIn}{A character vector specifying the motif sequences to be detected and quantified. Ambiguities
  can be specified using IUPAC codes or [ ] annotations. Optionally, g-quadruplexes can be specified as "G4".
  }
  \item{seqType}{A string specifying the type of sequence being provided to search for motifs. The options are,
  \code{"dna"}, \code{"rna"}, or \code{"protein"}. The default is \code{"dna"}.
  }
  \item{dist}{A numeric value specifying the minimal nucleotide distance between motifs. The default is 1.
  }
  \item{min_score}{A numeric value specifying the threshold for the quality of the G-quadruplexes prediction. 
  This is only required when \code{motifsIn} includes the "G4" option. The default is 47. See \link[pqsfinder]{pqsfinder}
  for details.
  }
  \item{unitOut}{A string to specify whether the output should be the \code{"number"} of motifs detected per gene
  (for a given sequence region), or the \code{"position"} of the detected motifs in the mRNA sequence (start and end). Note
  that only the output with \code{"number"} can be used in the downstream analysis with the \link[postNet]{featureIntegration} 
  function. The default is \code{"number"}.
  }
  \item{resid}{Logical indicating if the quantification of motifs should be corrected for the length of the
  sequences. If \code{TRUE}, the values returned are the residuals from the linear model: # of motifs ~ log2(length).
  Note that this option can only be \code{TRUE} when \code{unitOut = "number"}. The default is \code{FALSE}.
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
  performed. Note that plots will only be generated if \code{unitOut = "number"}. The default is \code{TRUE}. 
  }
}

\details{The number of motifs in the selected region for each gene is calculated. There is also an option to predict the presence of g-quadruplexes using the pqsfinder package https://bioconductor.org/packages/release/bioc/html/pqsfinder.html . If plotOut is TRUE, all indicated comparisons between different sets of genes are performed and plotted.
}

\value{List of vectors for selected regions with calculated number of motifs present for each gene. The function will also output plots with desired comparisons.
}

\references{
If you use the in-built "G4" option in your analysis, please cite pqsfinder: 

Hon J, Martínek T, Zendulka J, Lexa M. pqsfinder: an exhaustive and imperfection-tolerant search tool for potential quadruplex-forming sequences in R. Bioinformatics. 2017 Nov 1;33(21):3373-3379. doi: 10.1093/bioinformatics/btx413. PMID: 29077807.
}

\seealso{
\link[pqsfinder]{pqsfinder}
}

\examples{

}
