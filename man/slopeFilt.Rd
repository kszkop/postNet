\name{slopeFilt}
\alias{slopeFilt}
\title{Filtering \pkg{anota2seq} analysis results according to slopes for use in GSEA, GAGE, and GO term analyses
}

\description{The \code{slopeFilt} function identifies the genes with unrealistic models of changes in translation efficiency according to the slopes fitted in the per-gene \pkg{anota2seq} analysis of partial variance (APV) models. These genes are then filtered out of the analysis by providing the output of \code{slopeFilt} to the \code{genesSlopeFiltOut} argument of the \code{\link[postNet]{gseaAnalysis}}, \code{\link[postNet]{gseaPlot}}, \code{\link[postNet]{gageAnalysis}}, and \code{\link[postNet]{goAnalysis}} functions. See the \pkg{anota2seq} vignette for additional details on slope filtering.      
}

\usage{
slopeFilt(ads,
          regulationGen,
          contrastSel,
          minSlope = NULL,
          maxSlope = NULL)
}

\arguments{
  \item{ads}{An S4 object of class \code{Anota2seqDataSet} (resulting from the \code{\link[anota2seq]{anota2seqRun}}
  function in the \pkg{anota2seq} package). 
}
  \item{regulationGen}{The regulation effect measurement that will be taken from \code{ads}. This parameter can be
  either "translation", or "buffering". Note that slope filtering is not relevant for changes in mRNA abundance. 
}
  \item{contrastSel}{The contrast in \code{ads} where the regulation effect measurement corresponding to 
  \code{regulationGen} should be extracted from. This must be a single numeric value.
}
  \item{minSlope}{The minimum threshold to filter gene identifiers whose calculated slopes in anota2seq APV models are
  too small to be realistic. If \code{regulationGen = "translation"}, the recommended threshold is -1 (i.e., excludes
  gene identifiers with a slope < (-1)). If \code{regulationGen = "buffering"}, the recommended threshold is -2 (i.e.,
  excludes gene identifiers with a slope < (-2)). Default is \code{NULL}, which will automatically apply the recommended
  filtering threshold for the selected \code{regulationGen}.
}
  \item{maxSlope}{The maximum threshold to filter gene identifiers whose calculated slopes in anota2seq APV models are
  too large to be realistic. If \code{regulationGen = "translation"}, the recommended threshold is 2 (i.e., excludes
  gene identifiers with a slope > 2). If \code{regulationGen = "buffering"}, the recommended threshold is 1 (i.e.,
  excludes gene identifiers with a slope > 1). Default is \code{NULL}, which will automatically apply the recommended
  filtering threshold for the selected \code{regulationGen}. 
}
}

\details{When performing GSEA, GAGE, or GO term analysis using the output of an \pkg{anota2seq} analysis, it is often necessary to filter the input genes and log2 fold changes for the "translation" and "buffering" regulatory modes prior to performing enrichment analyses. This is due to that the slopes fitted by the anota2seq APV models can sometimes have unrealistic values, or suggest unlikely translational regulation, impacting the analysis of changes in translation or translational buffering (or offsetting). Filtering out genes with these unrealistic slopes is especially important for GSEA and GAGE analyses, which rely on rankings. For analyses relying on hypergeometric tests, such as GO term enrichment, the impact of filtering on the analysis is likely to be more negligible. However, slope filtering is still recommended. Note that in high-quality data sets, usually few genes will require slope filtering.   
}

\value{A character vector of the gene identifiers that will be excluded from downstream enrichment analyses based on the selected slope filtering thresholds.
}

\seealso{
\code{\link[anota2seq]{anota2seqSelSigGenes}} \cr
\code{\link[anota2seq]{anota2seqRun}} \cr
\code{\link[postNet]{gseaAnalysis}} \cr
\code{\link[postNet]{gseaPlot}} \cr
\code{\link[postNet]{gageAnalysis}} \cr
\code{\link[postNet]{goAnalysis}} \cr
\code{\link[postNet]{miRNAAnalysis}}
}

\examples{

# load example data:
 data("postNetExample", package = "postNetParcel")

# Initialize Anota2seqDataSet (see anota2seq vignette for details)
 ads <- anota2seq::anota2seqDataSetFromMatrix(
    dataP = postNetExample$ads_data$dataP,
    dataT = postNetExample$ads_data$dataT,
    phenoVec = postNetExample$ads_data$phenoVec,
    batchVec = c(1, 2, 3, 4, 1, 2, 3, 4),
    dataType = "RNAseq",
    normalize = FALSE)
    
# Run an anota2seq analysis:
# Note that the quality control and residual outlier testing are not 
# performed to limit the running time of this example. For full details 
# on running an analysis please see the anota2seq vignette and help manual.
 ads <- anota2seq::anota2seqRun(ads,
    performQC = FALSE, 
    performROT = FALSE, 
    useProgBar = FALSE)

# Get the genes to be filtered out of downstream enrichment analyses 
# using the buffering regulatory mode:
 filtOutGenes <- slopeFilt(ads,
                  regulationGen = "buffering",
                  contrastSel = 1)
                  
  str(filtOutGenes)
    
}
