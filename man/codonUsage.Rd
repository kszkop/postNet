\name{codonUsage}
\alias{codonUsage}
\title{Codon usage analysis
}

\description{The \code{codonUsage} function performs codon (or alternatively amino acid [AA]) usage analysis.
The analysis identifies single (or multiple, e.g. dicodons) codons enriched or depleted between different gene sets.
In the first step, for each gene set comparison, a chi-square test is applied to assess significant differences in
codon usage. Then, for significant codons, the odds ratio for each desired comparison pair is calculated and plotted 
against the codon frequency. The codons (or AA) of interest are identified as those with a high/low oddsratio and 
high frequency.
}

\usage{
codonUsage(ptn, 
          annotType = NULL, 
          sourceSeq = "load", 
          analysis, 
          codonN = 1, 
          dpn_path = NULL, 
          cds_filt = TRUE, 
          pAdj = 0.01, 
          rem5 = TRUE, 
          plotHeatmap = TRUE, 
          thresX1 = 0.3, 
          thresY1 = 0.3, 
          thresX2 = 0.3, 
          thresY2 = 0.3, 
          subregion = NULL, 
          subregionSel = NULL, 
          inputTable = NULL, 
          comparisons, 
          plotType_index = "boxplot", 
          pdfName = NULL)
}

\arguments{
  \item{ptn}{
  }
  \item{annotType}{select source of reference coding region sequences. Choose from "ccds", "custom", 
  "refseq". The default is "ccds". 
  }
  \item{sourceSeq}{if annotType is "ccds", choose "create" to automatically prepare reference sequences 
  from ccds database,  "load" to load preprepared sequences for selected species, or "custom" to use 
  customFile. The default is "load". 
  }
  \item{analysis}{select type of the analysis. Choose from "codon" or alternatively amino acid level 
  "AA". The default is "codon". 
  }
  \item{codonN}{n number of codon combinations. The default is 1. For di-codon choose 2 and so on. 
  The option is available only for type=="sequence". 
  }
  \item{dpn_path}{path to the dpn files obtained using the riboUtils function. It is required only 
  if type = "riboseq". 
  }
  \item{cds_filt}{
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
  \item{inputTable}{Dataframe containing names of dpn files in the first column, condition in numeric format 
  in the second column.
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
A table containing for each gene: the gene list it belongs to, as well as several codon indexes such as T3s (frequency of T at 3rd codon position), C3s (frequency of C at 3rd codon position), A3s, G3s, CAI (codon adaptation index), CBI (Codon Bias Index), Fop (Frequency of optimal codons), Nc (effective number of codons), GC3s (frequency of G or C at 3rd codon position), GC (GC content), L_sym (Frequency of synonymous codons), L_aa (Number of amino acids in protein), Gravy (general average hydropathicity score), Aromo (frequency of aromatic amino acids (Phe, Tyr, Trp) in the hypothetical translated gene product), tai (tRNA adaptation index)}}
 Sharp, P. M., and W. H. Li, (1987). The codon adaptation index a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research 15: 1281-1295.
 Bennetzen, J. L., and B. D. Hall, (1982). Codon selection in yeast. Journal of Biological Chemistry 257: 3026-3031.
 Ikemura, T., (1981). Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes: a proposal for a synonymous codon choice that is optimal for the E. coli system. Journal of Molecular Biology 151: 389-409.
 Wright, F., (1990). The effective number of codons used in a gene. Gene  87: 23-29.
 Kyte, J., and R. Doolittle, (1982). A simple method for displaying the hydropathic character of a protein. Journal of Molecular Biology 157: 105-132.
 Lobry, J. R., and C. Gautier, (1994). Hydrophobicity, expressivity and aromaticity are the major trends of amino acid usage in 999 Escherichia coli chromosome encoded genes. Nucleic Acids Research 22:3174-3180.
 dos Reis M, Savva R, Wernisch L. Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res. 2004 Sep 24;32(17):5036-44.
 dos Reis M, Wernisch L, Savva R. Unexpected correlations between gene expression and codon usage bias from microarray data for the whole Escherichia coli K-12 genome. Nucleic Acids Res. 2003 Dec 1;31(23):6976-85.
}

\value{List with selected codons in up and down sets for each desired comparison. In addition, summary codons calculation (list entry: codonAll) is also provided. The function will also output plots for desired comparisons.
}

\references{
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
