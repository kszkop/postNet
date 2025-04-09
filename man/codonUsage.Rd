\name{codonUsage}
\alias{codonUsage}
\title{Codon usage analysis
}

\description{The \code{codonUsage} function performs codon (or alternatively amino acid [AA]) usage analysis.
The analysis identifies single (or multiple, e.g. dicodons) codons enriched or depleted between different gene 
sets of interest. In the first step, for each gene set comparison, a chi-square test is applied to assess 
significant differences in usage. Then, for significant codons (or AA), the odds ratio for each desired 
comparison pair is calculated and plotted against the frequency. The codons (or AA) of interest are identified as 
those with a high/low odds ratio and high frequency.
}

\usage{
codonUsage(ptn, 
          annotType = NULL, 
          sourceSeq = "load", 
          analysis, 
          codonN = 1, 
          pAdj = 0.01, 
          rem5 = TRUE, 
          plotHeatmap = TRUE, 
          thresOddsUp = 0.3, 
          thresFreqUp = 0.3, 
          thresOddsDown = 0.3, 
          thresFreqDown = 0.3, 
          subregion = NULL, 
          subregionSel = NULL, 
          comparisons, 
          plotType_index = "boxplot", 
          pdfName = NULL)
}

\arguments{
  \item{ptn}{A \code{postNetData} object.
  }
  \item{annotType}{A string specifying the reference coding sequences to be used for the analysis. The options 
  are: \code{"ccds"} (the NCBI Consensus CDS database), or \code{"ptnCDS"} (the CDS sequences already stored in 
  the \code{postNetData} object. The default is \code{"ccds"}. 
  }
  \item{sourceSeq}{A string specifying the source of \code{annotType}. If \code{annotType = "ccds"} choose "create" 
  to automatically prepare reference sequences from ccds database,  "load" to load preprepared sequences for selected 
  species, or "custom" to use customFile. The default is "load". 
  }
  \item{analysis}{A string specifying thes type of the analysis to perform. The options are \code{"codon"}, or \code{"AA"} 
  to assess amino acid usage. The default is \code{"codon"}. 
  }
  \item{codonN}{An integer specifying the number of codon combinations to count. For example, \code{codonN = 1} will count each
  individual codon, while \code{codonN = 2}  will count di-codons, and so on. Note that this option is only available only 
  for \code{analysis = "codon"}. Note that when \code{codonN} > 1 additional indexes and plots will not be generated.
  The default is \code{1}. 
  }
  \item{pAdj}{A numeric value providing the adjusted p-value threshold for chi-square test results. The default is \code{0.01}. 
  }
  \item{rem5}{Logical specifying whether codons with fewer than five counts should be filtered out from the contingency table 
  (counts of codons by geneset) before performing the chi-square test. The default is \code{TRUE}. 
  }
  \item{plotHeatmap}{Logical specifying whether to plot a heatmap of the standardized residuals for each codon from the
  Chi-square test. The default is \code{TRUE}. 
  }
  \item{thresOddsUp}{select the top percentage threshold for odds ratio for selecting codons in up set. 
  The default 0.3 indicates that 30\% of codons with highest odds ratio is selected. It is combined 
  the thresFreqUp to select the most frequent codons with the highest odds ratio.
  }
  \item{thresFreqUp}{select the top percentage threshold for frequency for selecting codons in the up set. 
  The default 0.3 indicates that 30\% of the most frequent codons are selected. It is combined the 
  thresOddsUp to select the most frequent codons with the highest odds ratio. 
  }
  \item{thresOddsDown}{select the top percentage threshold for odds ratio for selecting codons in down set. 
  The default 0.3 indicates that 30\% of codons with lowest odds ratio is selected. It is combined 
  the thresFreqDown to select the most frequent codons with the lowest odds ratio.
  }
  \item{thresFreqDown}{select the top percentage threshold for frequency for selecting codons in the down set.
  The default 0.3 indicates that 30\% of the most frequent codons are selected. It is combined the 
  thresOddsDown to select the most frequent codons with the lowest odds ratio.
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
  \item{plotType_index}{
  }
  \item{pdfName}{Name to be appended to output PDF files. Default is \code{NULL}.
  }
}

\details{
}

\value{List with selected codons in up and down sets for each desired comparison. In addition, summary codons 
calculation (list entry: codonAll) is also provided. The function will also output plots for desired comparisons.
}

\references{
The NCBI Consensus CDS database:
\url{https://www-ncbi-nlm-nih-gov.proxy.kib.ki.se/projects/CCDS/CcdsBrowse.cgi}
}

\seealso{
}

\examples{
# If only the anota2seq object
codonOut <- codonUsage(annotType = "ccds", sourceSeq = "load", species = "human", analysis = "codon", codonN = 1, type = "sequence", pAdj = 0.01 , plotHeatmap = TRUE, ads = anota2seqOutput, regulation = c("translationUp","translationDown"), contrast = contrast = c(1,1), selection = "random", comparisons = list(c(0,1),c(0,2),c(1,2)))

# If only the geneList
codonOut <- codonUsage(annotType = "ccds", sourceSeq = "load", species = "human", analysis = "codon", codonN = 1, type = "sequence", pAdj = 0.01 , plotHeatmap = TRUE, geneList = inList,  selection = "random", comparisons = list(c(0,1),c(0,2),c(1,2)))

### For ribosomal profiling
#It is important that the annotation is the same as the one used for alignement. Optionally, it can be provided if different to the one in the package
annot <- retrieveFormatData(source='load', species='human')

dpn_files <- c(cond1_rep1.dpn, cond1_rep2.dpn, cond1_rep3.dpn,cond2_rep1.dpn, cond2_rep2.dpn, cond2_rep3.dpn)
condition <- c(1,1,1,2,2,2)
inputTable <- data.frame(dpnFiles=dpn_files, condition=condition,stringsAsFactors=FALSE)

codonOut <- codonUsage(type = "riboseq", annot=annot, dpn_path = '.', analysis = "codon", codonN = 1, pAdj = 0.01 , plotHeatmap = TRUE, ads = anota2seqOutput, regulation = c("translationUp","translationDown"), contrast = contrast = c(1,1), selection = "random", comparisons = list(c(0,1),c(0,2),c(1,2)), inputTable = inputTable)
}
