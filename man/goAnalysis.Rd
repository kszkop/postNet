\name{goAnalysis}
\alias{goAnalysis}
\title{Perform Gene Ontology Analysis on gene lists of interest with a \code{postNetData} object
}

\description{The \code{goAnalysis} function implements \link[clusterProfiler]{clusterProfiler} to perform Gene Ontology 
enrichment analysis using the gene sets associated with the regulatory effect measurement contained in a 
\code{postNetData} object. These gene sets are specified by the \code{regulationGen}, or \code{effectMeasure} parameters
in the \link[postNet]{postNetStart} function. 
}

\usage{
goAnalysis(ptn, 
          genesSlopeFiltOut = NULL, 
          category, 
          maxSize = 500, 
          minSize = 10, 
          counts = 10, 
          FDR = 0.15, 
          name = NULL)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
  \item{genesSlopeFiltOut}{If using an \code{Anota2seqDataSet} with analysis results for the "translation" or
  "buffering" regulatory modes, optionally provide a character vector with gene IDs to be excluded from GSEA that 
  is the output of the \code{\link[postNet]{slopeFilt}} function. Although optional, this step is strongly 
  recommended. This filtering is not necessary if using "mRNAAbundance", "totalmRNA", or "translatedmRNA". 
  Default is \code{NULL}.  
  }
  \item{category}{A character vector specifying one or more Gene Ontology categories to be included in the analysis. 
  The options are \code{'BP'}, \code{'CC'}, \code{'MF'}, and \code{'KEGG'}. For example, \code{category = c('BP','MF')}.
  }
  \item{maxSize}{Integer specifying the maximal size of a gene set to test. Terms/pathways with more genes will 
  be exluded. Default is \code{500}.
  }
  \item{minSize}{Integer specifying the minimal size of a gene set to test. Terms/pathways with fewer genes will 
  be exluded. Default is \code{10}.
  }
  \item{counts}{Integer specifying the minimal number of genes in the gene set of interest that must be included in the term.
  Terms with fewer genes overlapping with the gene set of interest than the threshold specified will be discarded, 
  regardless of significance. Default is \code{10}.
  }
  \item{FDR}{A numeric value providing the FDR threshold to be used to select significant enrichments. Default is
  \code{0.15}
  }
  \item{name}{Optionally, provide a name for the output file storing the results of the GO analysis. Note 
  that the "GO" slot of the \code{postNetData} object will reflect the results of the most recent
  analysis. However, multiple analyses can be run and saved by providing unique names to the output file. Default is
  \code{NULL}.
  }
}

\details{After running GO term analysis, the results can be retrieved from the \code{postNetData} object
using the S4 method \link[postNet]{ptn_GO} function. Results can also be visualized using the
\link[postNet]{goDotPlot} function.

Note that if the results of an \pkg{anota2seq} analysis are being used for GO term analysis, it is
recommended to first perform filtering using the \link[postNet]{slopeFilt} function.
}
\value{An S4 \link[DOSE]{enrichResult-class} object is stored in the \code{postNetData} object analysis slot for each category of GO
terms assessed. In addition, the "results" slot of each \code{enrichResult} object is written to an excel file. 
}
\references{
See \link[clusterProfiler]{clusterProfiler} for full details of enrichment analysis implentation in R. \cr

If you use the \code{goAnalysis} function in your analysis, please cite: \cr

S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang, W Xie, T Wu, L Xie, G Yu.
  Using clusterProfiler to characterize multiomics data. Nature Protocols. 2024, 19(11):3292-3320
}

\note{
}

\seealso{
\link[clusterProfiler]{clusterProfiler} \cr
\link[postNet]{postNetStart} \cr
\link[postNet]{slopeFilt} \cr
\link[postNet]{ptn_GO} \cr
\link[postNet]{goDotplot} \cr

}
\examples{
# ------------------------------------------------------------------
# Example 1: Running GO term enrichment analysis
# ------------------------------------------------------------------

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

# Run GO term analysis:
 ptn <- goAnalysis(ptn = ptn, 
                  category=c('BP'),
                  name = 'example')

# Extract the signficant enrichment results from the postNetData object:
 goOut <- ptn_GO(ptn, category = 'BP', geneList = 'translationUp', threshold = 0.05)
 str(goOut)
 
# Plot GO term analysis results:
 goDotplot(ptn=ptn, 
          category = 'BP', 
          nCategories=50, 
          pool=T, 
          size = 'geneRatio', 
          pdfName='example')
          
\dontrun{

# ------------------------------------------------------------------
# Example 2: Running GO term enrichment analysis with anota2seq results
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
 ads <- anota2seqRun(ads,
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
                          regulationGen="translation", 
                          contrastSel=1)
                    
# Create example custom gene sets for GSEA:
 set1 <- sample(row.names(postNetExample$ads_data$dataP), 15)
 inSet <- list(Set1 = set1)                   
                          
# Run GSEA on the C5 collection:
 ptn <- gseaAnalysis(ptn = ptn,
                    genesSlopeFiltOut = filtOutGenes,
                    geneSet = inSet,
                    name = 'ads_example')   
  }
}
