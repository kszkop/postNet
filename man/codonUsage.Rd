\name{codonUsage}
\alias{codonUsage}
\title{Codon usage analysis
}

\description{The \code{codonUsage} function performs codon (or alternatively amino acid [AA]) usage analysis.
The analysis identifies single (or multiple, e.g. dicodons) codons enriched or depleted between different gene 
sets of interest. In the first step, for each gene set comparison, a chi-square test is applied to assess 
significant differences in usage. Then, for significant codons, the odds ratio for each desired comparison 
pair is calculated and plotted against the codon frequency. The codons (or AA) of interest are identified as 
those with a high/low oddsratio and high frequency.
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
          thresX1 = 0.3, 
          thresY1 = 0.3, 
          thresX2 = 0.3, 
          thresY2 = 0.3, 
          subregion = NULL, 
          subregionSel = NULL, 
          comparisons, 
          plotType_index = "boxplot", 
          pdfName = NULL)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
  \item{annotType}{A string specifying the reference coding sequences to be used for the analysis. The options 
  are: \code{"ccds"} (the NCBI Consensus CDS database), or \code{"ptnCDS"} (the CDS sequences already stored in 
  the \code{postNetData} object. The default is \code{"ccds"}. 
  }
  \item{sourceSeq}{A string specifying the source of \code{annotType}. If \code{annotType = "ccds"} choose "create" to automatically prepare reference sequences 
  from ccds database,  "load" to load preprepared sequences for selected species, or "custom" to use 
  customFile. The default is "load". 
  }
  \item{analysis}{select type of the analysis. Choose from "codon" or alternatively amino acid level 
  "AA". The default is "codon". 
  }
  \item{codonN}{n number of codon combinations. The default is 1. For di-codon choose 2 and so on. 
  The option is available only for type=="sequence". 
  }
  \item{pAdj}{adjusted p-value threshold for chi-square test. The default is 0.01
  }
  \item{rem5}{
  }
  \item{plotHeatmap}{option to plot heatmap. The default is TRUE. 
  }
  \item{thresX1}{select the top percentage threshold for odds ratio for selecting codons in up set. 
  The default 0.3 indicates that 30\% of codons with highest odds ratio is selected. It is combined 
  the thresY1 to select the most frequent codons with the highest odds ratio.
  }
  \item{thresY1}{select the top percentage threshold for frequency for selecting codons in the up set. 
  The default 0.3 indicates that 30\% of the most frequent codons are selected. It is combined the 
  thresX1 to select the most frequent codons with the highest odds ratio. 
  }
  \item{thresX2}{select the top percentage threshold for odds ratio for selecting codons in down set. 
  The default 0.3 indicates that 30\% of codons with lowest odds ratio is selected. It is combined 
  the thresY2 to select the most frequent codons with the lowest odds ratio.
  }
  \item{thresY2}{select the top percentage threshold for frequency for selecting codons in the down set.
  The default 0.3 indicates that 30\% of the most frequent codons are selected. It is combined the 
  thresX2 to select the most frequent codons with the lowest odds ratio.
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
  \item{pdfName}{Name to be appended to output PDF files.
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
