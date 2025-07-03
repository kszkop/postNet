\name{codonUsage}
\alias{codonUsage}
\title{Analysis of differential codon usage
}

\description{The \code{codonUsage} function performs codon (or alternatively amino acid [AA]) usage analysis.
The analysis identifies single (or multiple, e.g. dicodons) codons enriched or depleted between different gene 
sets of interest. In addition, various codon usage indexes are compared between gene sets of interest, and differences in average counts and frequency are visualized as scatterplots and heatmaps. 
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
  \item{ptn}{A \linkS4class{postNetData} object.
  }
  \item{annotType}{A string specifying the reference coding sequences to be used for the analysis. The options 
  are: \code{"ccds"} (the NCBI Consensus CDS database), or \code{"ptnCDS"} (the CDS reference sequences already stored in 
  the \code{postNetData} object). The default is \code{"ccds"}. 
  }
  \item{sourceSeq}{A string specifying the source of \code{annotType}. This parameter is only required if \code{annotType = "ccds"}.
  The options are: \code{"create"} to automatically prepare reference sequences from ccds database, or \code{"load"} to load 
  pre-prepared sequences for selected species. The default is "load". 
  }
  \item{analysis}{A string specifying the type of the analysis to perform. The options are \code{"codon"}, or \code{"AA"} 
  to assess amino acid usage. The default is \code{"codon"}. 
  }
  \item{codonN}{An integer specifying the number of codon combinations to count. For example, \code{codonN = 1} will count
  each individual codon, while \code{codonN = 2}  will count di-codons, and so on. Note that this option is only available
  for \code{analysis = "codon"} and that when \code{codonN} > 1, odds ratios will not be calculated and additional indexes 
  and plots will not be generated. The default is \code{1}. 
  }
  \item{pAdj}{A numeric value specifying the adjusted p-value threshold for selecting significant Chi-square test results. 
  The default is \code{0.01}. 
  }
  \item{rem5}{Logical specifying whether codons with fewer than five counts should be filtered out from the contingency table 
  (counts of codons by gene set) before performing the Chi-square test. The default is \code{TRUE}. 
  }
  \item{plotHeatmap}{Logical specifying whether to plot a heatmap of the standardized residuals for each codon from the
  Chi-square test. The default is \code{TRUE}. 
  }
  \item{thresOddsUp}{A numeric value between 0 and 1 specifying the percentage threshold for the odds ratio in selecting
  codons or AA in the enriched set. For example, \code{thresOddsUp = 0.3} indicates that the top 30\% of enriched
  codons with the highest Odds ratios will be selected. The default is \code{0.3}.
  }
  \item{thresFreqUp}{A numeric value between 0 and 1 specifying the percentage threshold for the frequency in selecting
  codons or AA in the enriched set. For example, \code{thresFreqUp = 0.3} indicates that the top 30\% of
  the most frequently occurring enriched codons are selected. This metric is combined with \code{thresOddsUp} to select the
  most frequent codons with the highest Odds ratio. The default is \code{0.3}.  
  }
  \item{thresOddsDown}{A numeric value between 0 and 1 specifying the percentage threshold for the odds ratio in selecting
  codons or AA in the depleted set. For example, \code{thresOddsDown = 0.3} indicates that the top 30\% of
  depleted codons with the lowest Odds ratios will be selected. The default is \code{0.3}.
  }
  \item{thresFreqDown}{A numeric value between 0 and 1 specifying the percentage threshold for the frequency in selecting
  codons or AA in the depleted set. For example, \code{thresFreqDown = 0.3} indicates that the top 30\% of
  the most frequently occurring depleted codons are selected. This metric is combined with \code{thresOddsDown} to select
  the most frequent codons with the lowest Odds ratio. The default is \code{0.3}.
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
  \item{comparisons}{A list of numeric vectors indicating which gene sets defined in the \code{regulation} or \code{geneList} slots
  of the \code{postNetData} object should be used in pairwise comparisons. The numeric vectors should correspond to the indexes of
  the gene sets, where 0 denotes the background. For example, \code{list(c(0,2),c(1,2))} indicates that the second gene set in
  \code{regulation} or \code{geneList} will be compared to the background (0), and the first and second gene sets will be compared
  to each other. The default is \code{NULL}.
  }
  \item{plotType_index}{A string to indicate the type of output plot to generate for codon indexes. Options are,     \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"boxplot"}.
  }
  \item{pdfName}{Name to be appended to output PDF files. Default is \code{NULL}.
  }
}

\details{In the first step of the analysis, for each gene set comparison, a Chi-square test is applied to assess significant differences in codon or AA usage. Then, for significant codons or AA, the Odds ratio for each desired comparison pair is calculated and plotted against the frequency. The codons or AA of interest are usually identified as those with a high/low Odds ratio and high frequency. Codons or AAs of interest can then be further assess using the \link[postNet]{codonCalc} function (see example below). 

Several codon indexes are also compared between gene sets of interest, including the CAI (codon adaptation index), CBI (Codon Bias Index), FOP (frequency of optimal codons), L_aa (Number of amino acids in protein), and tAI (tRNA adaptation index).
}

\value{The results returned by the \code{codonUsage} function are compiled into an S4 object of class \code{postNetCodons} and stored in the "codons" slot of the \code{postNetData} object. The "codonAnalysis" slot provides a summary of calculations for all codons for each gene, and stores an S4 object of class \code{postNetCodonsAll} with the following slots: geneID, codon, AA, count, frequency, AACountPerGene, and relative_frequency. Here, frequency corresponds to the number of codons relative to all codons in the gene, and relative_frequency corresponds to the number of codons relative to all synonymous codons in the gene. The contents of the "codonAnalysis" slot can be accessed using the \link[postNet]{ptn_codonAnalysis} S4 method.     

The "codonSelection" slot stores a named list of the selected codons or AAs that were signficantly enriched or depleted (according to the selected thresholds) in the specified comparisons. This list can be used as input for the \link[postNet]{codonCalc} function (see example below). The contents of the "codonSelection" slot can be accessed using the \link[postNet]{ptn_codonSelection} S4 method. 

The \code{codonUsage} function also generates plots comparing the codon indexes described above between gene sets, a clustered heatmap of the Chi-square residuals for each codon between the gene sets of interest, and plots comparing the average codon count and frequency between the gene sets of interest (codons encoding the same amino acid are coloured and connected by lines).  
}

\references{
The NCBI Consensus CDS database: \url{https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi} \cr
\cr
CAI: \cr
Sharp, P. M., and W. H. Li, (1987). The codon adaptation index a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research 15: 1281-1295. \cr
\cr
CBI: \cr
Bennetzen, J. L., and B. D. Hall, (1982). Codon selection in yeast. Journal of Biological Chemistry 257: 3026-3031. \cr
\cr
FOP: \cr
Ikemura, T., (1981). Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes: a proposal for a synonymous codon choice that is optimal for the E. coli system. Journal of Molecular Biology 151: 389-409. \cr
\cr
tAI: \cr
dos Reis M, Savva R, Wernisch L. Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Res. 2004 Sep 24;32(17):5036-44. \cr
}

\seealso{
\link[postNet]{codonCalc} \cr
\link[postNet]{ptn_codonAnalysis} \cr
\link[postNet]{ptn_codonSelection} \cr
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
 
# Run codon usage analysis of single codons
 ptn <- codonUsage(ptn = ptn, 
          annotType = 'ptnCDS', 
          sourceSeq = "load", 
          analysis = 'codon', 
          codonN = 1, 
          pAdj = 0.01, 
          rem5 = TRUE, 
          plotHeatmap = FALSE, 
          thresOddsUp = 0.3, 
          thresFreqUp = 0.3, 
          thresOddsDown = 0.3, 
          thresFreqDown = 0.3, 
          subregion = NULL, 
          subregionSel = NULL, 
          comparisons = list(c(1,2)), 
          plotType_index = "violin", 
          pdfName = 'codons')
          
# Select codons of interest that were significantly enriched or depleted with high frequency, and the highest and lowest Odds ratios
 
 codons <- ptn_codonSelection(ptn, comparison = 1)
 
# Count and plot ecdfs for sets of enriched/depleted codons, and compare between gene sets of interest

 codonCounts  <- codonCalc(ptn = ptn,
                          analysis='codon', 
                          featsel= codons,
                          unit='count', 
                          comparisons = list(c(1,2)),
                          plotType = 'ecdf')
                      
 str(codonCounts)
}
