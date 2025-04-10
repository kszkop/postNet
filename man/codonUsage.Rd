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
  \item{thresOddsUp}{A numeric value between 0 and 1 specifying the percentage threshold for the odds ratio in selecting
  codons (or amino acids) in enriched set. For example, \code{thresOddsUp = 0.3} indicates that the top 30\% of enriched
  codons with the highest odds ratios will be selected. The default is \code{0.3}.
  }
  \item{thresFreqUp}{select the top percentage threshold for frequency for selecting codons in the up set. 
  The default 0.3 indicates that 30\% of the most frequent codons are selected. It is combined the 
  thresOddsUp to select the most frequent codons with the highest odds ratio. 
  }
  \item{thresOddsDown}{A numeric value between 0 and 1 specifying the percentage threshold for the odds ratio in selecting
  codons (or amino acids) in the depleted set. For example, \code{thresOddsUp = 0.3} indicates that the top 30\% of depleted
  codons with the highest odds ratios will be selected. The default is \code{0.3}.
  }
  \item{thresFreqDown}{select the top percentage threshold for frequency for selecting codons in the down set.
  The default 0.3 indicates that 30\% of the most frequent codons are selected. It is combined the 
  thresOddsDown to select the most frequent codons with the lowest odds ratio.
  }
  \item{subregion}{Optionally, it is possible to specify a more specific subregion of the sequences to either
  select or exclude from the analysis. This is done by providing a numeric value indicating the number of
  nucleotides from either the start of the sequence region if positive, or the end if negative. For example,
  \code{subregion = 50}, would denote the first 50 nucleotides of the CDS, while \code{subregion = -50}, would 
  denote the last 50 nucleotides of the CDS. The default is \code{NULL}.
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
  \item{plotType_index}{If \code{plotOut} is \code{TRUE}, provide a string to indicate the type of output 
  plot to generate for codon indexes. Options are, \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. 
  The default is \code{"boxplot"}.
  }
  \item{pdfName}{Name to be appended to output PDF files. Default is \code{NULL}.
  }
}

\details{A table containing for each gene: the gene list it belongs to, as well as several codon indexes such as CAI (codon adaptation index), CBI (Codon Bias Index),
Fop (Frequency of optimal codons), GC3s (frequency of G or C at 3rd codon position), L_aa (Number of amino acids in protein), tai (tRNA adaptation index)

}

\value{List with selected codons in up and down sets for each desired comparison. In addition, summary codons 
calculation (list entry: codonAll) is also provided. The function will also output plots for desired comparisons.
}

\references{
The NCBI Consensus CDS database:
\url{https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi}

Sharp, P. M., and W. H. Li, (1987). The codon adaptation index a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research 15: 1281-1295.

Bennetzen, J. L., and B. D. Hall, (1982). Codon selection in yeast. Journal of Biological Chemistry 257: 3026-3031.

Ikemura, T., (1981). Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes: a proposal for a synonymous codon choice that is optimal for the E. coli system. Journal of Molecular Biology 151: 389-409.

Bernardi G, Olofsson B, Filipski J, Zerial M, Salinas J, Cuny G, Meunier-Rotival M, Rodier F. The mosaic genome of warm-blooded vertebrates. Science. 1985 May 24;228(4702):953-8. doi: 10.1126/science.4001930. PMID: 4001930.

Wright, F., (1990). The effective number of codons used in a gene. Gene  87: 23-29.

Kyte, J., and R. Doolittle, (1982). A simple method for displaying the hydropathic character of a protein. Journal of Molecular Biology 157: 105-132.

Lobry, J. R., and C. Gautier, (1994). Hydrophobicity, expressivity and aromaticity are the major trends of amino acid usage in 999 Escherichia coli chromosome encoded genes. Nucleic Acids Research 22:3174-3180.

dos Reis M, Savva R, Wernisch L. Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res. 2004 Sep 24;32(17):5036-44.

dos Reis M, Wernisch L, Savva R. Unexpected correlations between gene expression and codon usage bias from microarray data for the whole Escherichia coli K-12 genome. Nucleic Acids Res. 2003 Dec 1;31(23):6976-85.
}

\seealso{
\link[postNet]{codonCalc} \cr
\link[postNet]{ptn_codonAnalysis} \cr
\link[postNet]{ptn_codonSelection} \cr
}

\examples{
}
