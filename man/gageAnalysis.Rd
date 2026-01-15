\name{gageAnalysis}
\alias{gageAnalysis}
\title{Perform GAGE analysis with a \code{postNetData} object 
}
\description{The \code{gageAnalysis} function implements the \link[gage]{gage} package to perform Generally Applicable Gene-set Enrichment (GAGE) analysis using the regulatory effect measurement contained in a \code{postNetData} object. These values are determined by the \code{regulationGen} or \code{effectMeasure} parameters in the \link[postNet]{postNetStart} function. 
}

\usage{
gageAnalysis(ptn,
            category, 
            genesSlopeFiltOut = NULL, 
            maxSize = 500, 
            minSize = 10)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }
  \item{category}{A character vector specifying one or more gene ontology categories to be included in the analysis. The options are \code{"BP"}, \code{"CC"}, \code{"MF"}, and \code{"KEGG"}. For example, \code{category = c("BP", "MF")}.
  }
  \item{genesSlopeFiltOut}{If using an \code{Anota2seqDataSet} with analysis results for the "translation" or "buffering" regulatory modes, optionally provide a character vector with gene IDs to be excluded from GAGE analysis. This can be generated using the \link[postNet]{slopeFilt} function. Although optional, this step is strongly recommended. This filtering is not necessary if using "mRNAAbundance", "totalmRNA", or "translatedmRNA". The default is \code{NULL}.
  }
  \item{maxSize}{Integer specifying the maximal size of a gene set to test. Terms/pathways with more genes will be excluded. The default is \code{500}.
  }
  \item{minSize}{Integer specifying the minimal size of a gene set to test. Terms/pathways with fewer genes will be excluded. The default is \code{10}.
  }
}

\details{After running GAGE analysis, the results can be retrieved from the \code{postNetData} object using the S4 method \link[postNet]{ptn_GAGE}.

Note that if using results from an \pkg{anota2seq} analysis, it is strongly recommended to first apply filtering using the \link[postNet]{slopeFilt} function to exclude genes with unrealistic regression slopes.
}

\value{The results returned by the \code{gageAnalysis} function are compiled into an S4 object of class \code{postNetGAGE} and stored in the "GAGE" slot of the \code{postNetData} object. As described in the \link[gage]{gage} package, each gene ontology category has a slot in the results object that stores a named list with three elements ("greater", "less" and "stats"). The "stats" element contains the test statistics. The "greater" and "less" elements are data.frames with the results of a two-directional test, each having the following columns:

 \itemize{ 
  \item \strong{p.geomean}: The geometric mean of the individual p-values from multiple single array-based gene set tests.
  
  \item \strong{stat.mean}: The mean of the individual statistics from multiple single array-based gene set tests. The value denotes the magnitude of the gene-set level changes, and the sign denotes the direction of the changes.
  
  \item \strong{p.val}: The global p-value or summary of the individual p-values from multiple single array-based gene set tests.
  
  \item \strong{q.val}: The Benjamini–Hochberg adjusted p-value, as implemented by the \pkg{multtest} package.
  
  \item \strong{set.size}: The number of genes included in the gene set.  
  
  \item \strong{Genes}: The gene IDs in the input data belonging to the gene ontology gene set.
  }
}

\references{
If you use the \code{gageAnalysis} function in your analysis, please cite: \cr

Luo, W., Friedman, M., Shedden K., Hankenson K., and Woolf P. GAGE: Generally Applicable Gene Set Enrichment for Pathways Analysis. BMC Bioinformatics 2009, 10:161.
}

\seealso{
\link[gage]{gage} \cr
\link[postNet]{ptn_GAGE} \cr 
\link[postNet]{slopeFilt} 
}

\examples{
# ------------------------------------------------------------------
# Example 1: Running GAGE with gene lists 
# ------------------------------------------------------------------

# load example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Run GAGE:
 ptn <- gageAnalysis(ptn,
                    category = "CC")

# Extract the significant enrichment results from the postNetData object:
 gageOut <- ptn_GAGE(ptn = ptn,
                    category = "CC",
                    direction = "greater", 
                    threshold = 1) # This will return the complete results. 
                                   # For only significant results 
                                   # the threshold should be lowered. 
 str(gageOut)

\dontrun{

# ------------------------------------------------------------------
# Example 2: Running GAGE requiring anota2seq analysis slope filtering
# ------------------------------------------------------------------

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
    
# Initialize the postNetData object, where "ads" is an Anota2seqDataSet object
# resulting from an anota2seqRun call:
  ptn <- postNetStart(
      ads = ads,
      regulation = c("translationUp","translationDown"),
      contrast = c(1,1),
      regulationGen = "translation",
      contrastSel = 1,
      selection = "random",
      setSeed = 123,
      source = "load",
      species = "human"
    )
    
# Run slope filtering to remove the genes with unrealistic slopes:
 filtOutGenes <- slopeFilt(ads = ads, 
                          regulationGen = "translation", 
                          contrastSel = 1)
                          
# Run GAGE:
 ptn <- gageAnalysis(ptn,
                    genesSlopeFiltOut = filtOutGenes,
                    category = "CC")   
  }
}
