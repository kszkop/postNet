\docType{methods}
\name{ptn_miRNA_to_gene}
\alias{ptn_miRNA_to_gene,postNetData-method}
\title{Accessor for the "miRNA_to_gene" slot of a \code{postNetData} object
}

\description{Retrieve the "miRNA_to_gene" slot from a \code{postNetData} object, which stores a list of genes that are predicted to be targeted by the miRNA that passed 
the filtering thresholds for \code{contextScore} and/or \code{Pct}. 
}

\usage{
\S4method{ptn_miRNA_to_gene}{postNetData}(ptn, miRNAs)
}

\arguments{
  \item{ptn}{A \code{postNetData} object.
  }
  \item{miRNAs}{A character vector specifying one or more miRNAs for which to extract the list of target genes. 
  }
}

\value{A list where each element is a vector of genes that are predicted to be targeted by the miRNA that passed the filtering thresholds for \code{contextScore} and/or \code{Pct}. 
}

\references{
If you use the \code{miRNAanalysis} function in your analysis, please cite:

McGeary SE, Lin KS, Shi CY, Pham T, Bisaria N, Kelley GM, Bartel DP. The biochemical basis of microRNA targeting efficacy. Science Dec 5, (2019). 

Agarwal V, Bell GW, Nam J, Bartel DP. Predicting effective microRNA target sites in mammalian mRNAs. eLife, 4:e05005, (2015). eLife Lens view.

Luo, W., Friedman, M., Shedden K., Hankenson, K. and Woolf, P GAGE: Generally Applicable Gene Set Enrichment for Pathways Analysis. BMC Bioinformatics 2009, 10:161.
}

\seealso{
\link[gage]{gage} \cr
\link[postNet]{miRNAanalysis} \cr
\link[postNet]{ptn_miRNA_analysis} \cr
\link[postNet]{postNetStart} \cr
}

\examples{
# load example data:
 postNetExample <- postNetParcel::postNetExample
 ptn <- postNetExample$ptn

# Run GAGE miRNA target enrichment analysis:
 ptn <- miRNAanalysis(ptn = ptn,
                     miRNATargetScanFile = 'miRNA_predictions_hsa.txt',
                     contextScore = -0.2,
                     Pct = 0.4)

# Extract the significant enrichment results from the postNetData object:
 miRNAtargets <- ptn_miRNA_to_gene(ptn = ptn,
                    miRNAs = c('hsa-miR-138-5p','hsa-miR-182-5p'))
 str(miRNAtargets)
}

