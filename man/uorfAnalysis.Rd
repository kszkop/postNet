\name{uorfAnalysis}
\alias{uorfAnalysis}
\title{Detect and quantify upstream open reading frames (uORFs)
}

\description{The \code{uorfAnalysis} function detects the presence and position of uORFs in 5'UTRs. Options are available to examine uORFs with both cannonical and non-cannonical start codons, as well as different Kozak contexts. In addition, the proportions of mRNA transcripts with uORFs can be plotted, and statistical comparisons can be made between gene sets of interest against background, and/or other genes sets.
}

\usage{
uorfAnalysis(ptn, 
            startCodon = "ATG", 
            KozakContext = "strong", 
            onlyUTR5 = FALSE, 
            unitOut = "number", 
            comparisons = NULL, 
            plotOut = TRUE, 
            pdfName = NULL)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
  \item{startCodon}{A string specifying the start codon. Only uORFs with the specified start codon will be detected.
  The default is the canonical start codon, \code{"ATG"}.
  }
  \item{KozakContext}{A string to select the Kozak Context. The options are: \code{"strong"} ([AG]..startCodonG),
  \code{"adequate1"} ([AG]..startCodon), \code{"adequate2"} (...startCodonG), \code{"weak"} (...startCodon.), or
  \code{"any"} (the nucleotide context of the start codon is disregarded). The default is \code{"strong"}. 
  }
  \item{onlyUTR5}{A logical to select whether the uORFs detected should be completely contained in the 5'UTR,
  (i.e., the stop codons of the uORFs occur before the start of the CDS). If \code{TRUE}, uORFs must be completely
  5'UTR-contained to be detected. If \code{FALSE}, only the start codon of the uORF must be in 5'UTR, with the
  possibility that the stop codon may be in the CDS or 3'UTR. The default is \code{TRUE}.
  }
  \item{unitOut}{A string to specify whether the output should be the \code{"number"} of uORFs detected per gene, or
  the \code{"position"} of the detected uORFs in the mRNA sequence (nucleotide position of the start and stop). Note
  that only the output with \code{"number"} can be used in the downstream analysis with the
  \link[postNet]{featureIntegration} function. The default is \code{"number"}.
  }
  \item{comparisons}{A list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the 
  values specified with \code{unitOut} will be returned, and statistical comparisons between gene sets will not be
  performed. Note that plots will only be generated if \code{unitOut = "number"}. The default is \code{TRUE}. 
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{A two-sided Wilcoxon Rank Sum test is performed to identify significant differences in uORFs detected between gene sets of interest, or against the background gene set.
}

\value{If the \code{unitOut = "number"}, the output will be a named list of vectors with the number of uORFs 
detected in a particular Kozak context for each gene. This list can be used with the downstream 
\link[postNet]{featureIntegration} function. If the \code{unitOut = "position"}, the output will be a named
list of lists with the start and end positions of each uORF for each gene. These positions cannot be used
in \link[postNet]{featureIntegration}. The \code{uorfAnalysis} function can also return PDF files of output
plots with statistical comparisons between gene sets of interest when \code{unitOut = "number"}.
}

\examples{
# load example data:
 postNetExample <- postNetParcel::postNetExample
 ptn <- postNetExample$ptn

# Identify uORFs with canonical start codons in a strong Kozak context, fully contained within 5'UTRs, 
# and compare between translationUp vs. translationDown genes

 uORFs <- uorfAnalysis(ptn = ptn,
                      comparisons = list(c(1,2)),
                      startCodon = "ATG",
                      KozakContext = c("strong"),
                      onlyUTR5 = TRUE,
                      unitOut = "number",
                      plotOut = TRUE,
                      pdfName = 'example')
                      
 str(uORFs)
}
