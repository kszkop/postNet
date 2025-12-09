\docType{function}
\name{contentMotifs}
\alias{contentMotifs}
\title{Detect, quantify, and compare sequence motifs between gene sets
}

\description{The \code{contentMotifs} function quantifies the number or position of motifs in sequence regions of interest and performs statistical comparisons between gene sets of interest, or against background. Motifs can be specified directly, or the output of the \link[postNet]{motifAnalysis} function, which detects de-novo motifs using STREME (part of the MEME-Suite), can be supplied. Finally, the function can also implement \link[pqsfinder]{pqsfinder} to identify G-quadruplexes. 
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
  \item{ptn}{A \linkS4class{postNetData} object.
  }
  \item{motifsIn}{A character vector specifying the motif sequences to be detected and quantified. Ambiguities can be specified using IUPAC codes or [ ] (bracket) annotations. Motifs discovered with \link[postNet]{motifAnalysis} can be retrieved using \link[postNet]{ptn_motifSelection} and used. Optionally, G-quadruplexes can be specified as \code{"G4"}.
  }
  \item{seqType}{A string specifying the type of sequence being provided to search for motifs. The options are \code{"dna"}, \code{"rna"}, or \code{"protein"}. The default is \code{"dna"}.
  }
  \item{dist}{A numeric value specifying the minimal distance between motifs (in nucleotides for DNA/RNA or amino acids for protein). The default is 1.
  }
  \item{min_score}{A numeric value specifying the threshold for the quality of the G-quadruplexes prediction. This is only required when \code{motifsIn} includes the "G4" option. The default is 47. See \link[pqsfinder]{pqsfinder} for details.
  }
  \item{unitOut}{A string to specify whether the output should be the \code{"number"} of motifs detected per gene (for a given sequence region), or the \code{"position"} of the detected motifs in the mRNA sequence (start and end). Note that only the output with \code{"number"} can be used in the downstream analysis with the \link[postNet]{featureIntegration} function. The default is \code{"number"}.
  }
  \item{resid}{Logical indicating if the quantification of motifs should be corrected for the length of the sequences. If \code{TRUE}, the values returned are the residuals from the linear model: number of motifs ~ log2(sequence region length). Note that this option can only be set to \code{TRUE} when \code{unitOut = "number"}. See details below for usage recommendations. The default is \code{FALSE}.
  }
  \item{region}{A character vector specifying the sequence region(s) to be analyzed. This can be \code{"UTR5"}, \code{"UTR3"}, \code{"CDS"}, or any combination of the regions.
  }
  \item{subregion}{Optionally, it is possible to specify a more specific subregion of the sequences to either select or exclude from the analysis. This is done by providing a numeric value indicating the number of nucleotides from either the start of the sequence region if positive, or the end if negative. For example, \code{subregion = 50} denotes the first 50 nucleotides of the sequence region(s) specified by \code{region}, while \code{subregion = -50} denotes the last 50 nucleotides of the sequence region(s). The default is \code{NULL}.
  }
  \item{subregionSel}{If a value is provided for \code{subregion}, indicate whether the subregion should be selected (limiting the analysis to this region) or excluded by specifying either \code{"select"} or \code{"exclude"}. The default is \code{NULL}.
  }
  \item{comparisons}{A list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}.
  }
  \item{pdfName}{Name to be appended to output PDF files. The default is \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the values specified with \code{unitOut} will be returned, and statistical comparisons between gene sets will not be performed. Note that plots will only be generated if \code{unitOut = "number"}. The default is \code{TRUE}. 
  }
}

\details{A two-sided Wilcoxon Rank Sum test is performed to identify significant differences in motif content between gene sets of interest, or against the background gene set. 

Note that the \code{contentMotifs} function uses a more strict method of sequence matching to count motifs than the MEME-Suite, and for this reason some motifs identified with \link[postNet]{motifAnalysis}, which implements STREME, may have divergent quantifications between methods. Using a more stringent threshold with the \code{stremeThreshold} parameter in \link[postNet]{motifAnalysis} may remedy this. Alternatively, motifs identified with STREME can be counted using the \link[memes]{runFimo} tool provided with the MEME-Suite.

If the output motifs will be used in downstream analysis with \link[postNet]{featureIntegration}, consider carefully whether to apply the \code{resid} correction. If the length of the sequence regions will be included in the \code{featureIntegration} model, it is not recommended to correct the motif content for the sequence length. However, if length will not be included in \code{featureIntegration}, the correction can be performed.
}

\value{If \code{unitOut = "number"}, the output will be a named list of vectors with the number of motifs detected for each gene in each sequence region specified. This list can be used as input for the downstream \link[postNet]{featureIntegration} analysis. If \code{unitOut = "position"}, the output will be a named list of lists with the start and end positions of each motif for each gene in each sequence region specified. These positions cannot be used in \link[postNet]{featureIntegration}. The \code{contentMotifs} function can also return PDF files of output plots with statistical comparisons between gene sets of interest when \code{unitOut = "number"}.
}

\references{
If you use the in-built "G4" option in your analysis, please cite pqsfinder: 

Hon J, Martínek T, Zendulka J, Lexa M. pqsfinder: an exhaustive and imperfection-tolerant search tool for potential quadruplex-forming sequences in R. Bioinformatics. 2017 Nov 1;33(21):3373-3379. doi: 10.1093/bioinformatics/btx413. PMID: 29077807.
}

\seealso{
\link[pqsfinder]{pqsfinder} \cr
\link[memes]{runFimo} \cr
\link[postNet]{motifAnalysis}
}

\examples{

# load example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn
 
# ------------------------------------------------------------------
# Example 1: Detect and quantify a given motif within 5'UTRs
# ------------------------------------------------------------------

# Detect and quantify a given motif within 5'UTRs, and compare between translationUp vs. translationDown genes

 UTR5_SGCSGCS_num <- contentMotifs(ptn = ptn, 
                                   motifsIn = "SGCSGCS",
                                   region = c("UTR5"), 
                                   comparisons = list(c(1,2)),
                                   dist = 1,
                                   unitOut = "number",
                                   pdfName = "example",
                                   plotOut = TRUE)
                      
 str(UTR5_SGCSGCS_num)
 
 # Now find the positions of the motif in the 5'UTR
 
 UTR5_SGCSGCS_pos <- contentMotifs(ptn = ptn, 
                                   motifsIn = "SGCSGCS",
                                   region = c("UTR5"), 
                                   comparisons = list(c(1,2)),
                                   dist = 1,
                                   unitOut = "position")
                      
 str(UTR5_SGCSGCS_pos)
 
 \dontrun{
# ------------------------------------------------------------------
# Example 2: Detect and quantify motifs identified using motifAnalysis
# ------------------------------------------------------------------

# First run motifAnalysis to identify significantly enriched motifs in the 3'UTR
 ptn <- motifAnalysis(ptn = ptn, 
                      stremeThreshold = 0.05,
                      minwidth = 6,
                      memePath = "/meme/bin", 
                      region = c("UTR3"))

# Quantify the presence of these motifs in transcripts, and compare between translationUp vs. translationDown genes, adjusting for the length of the sequence.
 
 UTR3_denovo_motifs <- contentMotifs(ptn = ptn, 
                                     motifsIn = ptn_motifSelection(ptn, region = "UTR3"),
                                     region = c("UTR3"), 
                                     comparisons = list(c(1,2)),
                                     dist = 1,
                                     resid = TRUE,
                                     unitOut = "number",
                                     pdfName = "example",
                                     plotOut = TRUE)
                      
 str(UTR3_denovo_motifs)
 
# ------------------------------------------------------------------
# Example 3: Detect and quantify potential G-quadruplex-forming sequences
# ------------------------------------------------------------------

# Quantify the presence of G4 motifs in 5'UTR and CDS, and compare between translationUp vs. translationDown 
# genes.
 
 G4_UT5_CDS <- contentMotifs(ptn = ptn, 
                             motifsIn = "G4",
                             region = c("UTR5","CDS"), 
                             comparisons = list(c(1,2)),
                             dist = 1,
                             resid = FALSE,
                             unitOut = "number",
                             pdfName = "example",
                             plotOut = TRUE)
                      
 str(G4_UT5_CDS)
 }
}
