\docType{methods}
\name{ptn_miRNA_analysis}
\alias{ptn_miRNA_analysis,postNetData-method}
\title{Accessor for the "miRNA_analysis" slot of a \code{postNetData} object
}

\description{Retrieve the "miRNA_analysis" slot from a \code{postNetData} object, which holds the results of GAGE analysis identifying miRNA target enrichments in
gene sets of interest created by the \link[postNet]{miRNAAnalysis} function.
}

\usage{
\S4method{ptn_miRNA_analysis}{postNetData}(ptn, 
                  direction, 
                  threshold)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
  \item{direction}{A character indicating the directionality of regulation of the gene sets. Can be either
  \code{"greater"}, or \code{"less"}.
  }
  \item{threshold}{A numeric value providing the FDR threshold to be used to select results from the GAGE miRNA enrichment analysis. 
  }
}

\details{Note that note that enrichments in miRNA predicted to target genes that are upregulated (e.g., if the regulatory effect measurement is log2 fold change) that appear
in the results table labelled "greater" can be interpretted as those miRNA that may be dwonregulated or otherwise not active in the experimental condition. Likewise, 
enrichments in miRNA predicted to target genes that are downregulated that appear in the results table labelled "less" can be interpretted as those miRNA that may be
upregulated or active in the experimental condition. 
}

\value{The value returned by the \code{ptn_miRNA_analysis} function is a data.frame with the following columns (see the \link[gage]{gage} package for full details):

 \itemize{ 
  \item \strong{p.geomean}: The geometric mean of the individual p-values from multiple single array-based gene set 
  tests.
  \item \strong{stat.mean}: The mean of the individual statistics from multple single array-based gene set tests. The
  value denotes the magnitude of the gene-set level changes, and the sign denotes the direction of the changes.
  \item \strong{p.val}: The global p-value or summary of the individual p-values from multiple single array-based gene
  set tests.
  \item \strong{q.val}: The BH adjusted p-value, as implemented by the \pkg{multtest} package.
  \item \strong{set.size}: The number of genes included in the gene set.
  }
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
\link[postNet]{ptn_miRNA_to_gene} \cr
\link[postNet]{postNetStart} \cr
}

\examples{
# load example data:
 data(postNetExample)
    
# Initialize the postNetData object with custom gene lists:
 set1 <- postNetExample$geneList[[1]][1:100]
 set2 <- postNetExample$geneList[[2]][1:100]
 myGenes <- list('Set1' = set1, 'Set2' = set2)
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

# Run GAGE miRNA target enrichment analysis:
 ptn <- miRNAanalysis(ptn,
                     miRNATargetScanFile = 'miRNA_predictions_hsa.txt',
                     contextScore = -0.2)

# Extract the signficant enrichment results from the postNetData object:
 miRNAout <- ptn_miRNA_analysis(ptn = ptn,
                    direction = "less", 
                    threshold = 1)
 str(miRNAout)
}
