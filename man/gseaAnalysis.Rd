\name{gseaAnalysis}
\alias{gseaAnalysis}
\title{Perform Gene Set Enrichment Analysis (GSEA) with a \code{postNetData} object 
}

\description{The \code{gseaAnalysis} function performs Gene Set Enrichment Analysis using the regulatory effect
measurement contained in a \code{postNetData} object. These values are specified by the \code{regulationGen}, or 
\code{effectMeasure} parameters in the \link[postNet]{postNetStart} function. 
}

\usage{
gseaAnalysis(ptn, 
             genesSlopeFiltOut = NULL, 
             collection = NULL, 
             subcollection = NULL, 
             subsetNames = NULL, 
             geneSet = NULL, 
             maxSize = 500, 
             minSize = 10, 
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
  \item{collection}{A character vector selecting the MSigDB collection of interest, for example:  \code{'c8'}, 
  or \code{'h'}. Although possible, it is not recommended to run all collections in the same analysis. Available
  options for collections can be found at \url{https://www.gsea-msigdb.org/gsea/msigdb/index.jsp}, or by using 
  the \pkg{msigdb} package (see example below). Currently, collections are only available for human and mouse. 
  }
  \item{subcollection}{Optionally, provide a character vector selecting the MSigDB subcollection of interest. 
  For example: \code{"GO:BP"}, or \code{"CP:REACTOME"}. Available options for subcollections can be found at 
  \url{https://www.gsea-msigdb.org/gsea/msigdb/index.jsp}, or by using the \pkg{msigdb} package (see example below).
  Default is \code{NULL}. 
  }
  \item{subsetNames}{Optionally, provide a character vector specifying specific terms to use for GSEA by providing 
  the names of the pathways. For example: \code{c("GOBP_RIBOSOME_BIOGENESIS"
  ,"GOBP_CAP_DEPENDENT_TRANSLATION_INITIATION")}. Default is \code{NULL}. 
  }
  \item{geneSet}{Optionally,instead using a GSEA collection, provide named list of character vectors with custom
  gene signatures to be used in GSEA. Default is \code{NULL}.
  }
  \item{maxSize}{Integer specifying the maximal size of a gene set to test. Terms/pathways with more genes will 
  be exluded. Default is \code{500}.
  }
  \item{minSize}{Integer specifying the minimal size of a gene set to test. Terms/pathways with fewer genes will 
  be exluded. Default is \code{10}.
  }
  \item{name}{Optionally, provide a name for the tab-delimited out file storing the results of the GSEA. Note 
  that the "GSEA" slot of the \code{postNetData} object will reflect the results of the most recent
  analysis. However, multiple analyses can be run and saved by providing unique names to the output file. Default is
  \code{NULL}.
  }
}

\details{After running GSEA, the results of the analysis can be retrieved from the \code{postNetData} object
using the S4 method \link[postNet]{ptn_GSEA} function. These results can then be provided as input to the
\link[postNet]{gseaPlot} function, which can be used to generate GSEA plots for all results, or for specific 
terms of interest. 
}

\value{A data.frame of GSEA results that is stored in the GSEA slot of the \code{postNetData} object, and written to a tab delimited text file. Each row of the data.frame corresponds to a term, with the following columns (also see \link[fgsea]{fgseaMultilevel} from the \pkg{fgsea} package):

 \itemize{ 
  \item \strong{Term}: The name of the term or pathway tested.
  
  \item \strong{ES}: The enrichment score.
  
  \item \strong{NES}: The normalized enrichment score, adjusted to the mean enrichment of randomly sampled 
  size-matched sets.
  
  \item \strong{log2err}: The expected error for the standard deviation of the P-value logarithm.
  
  \item \strong{count}: The number of genes in the set belonging to the term or pathway. 
  
  \item \strong{size}: The size of the pathway after removing genes not present in the input.  
  
  \item \strong{pvalue}: The enrichment p-value.  
  
  \item \strong{adjusted_pvalue}: The BH-adjusted p-value.
  
  \item \strong{Genes}: The gene IDs in the input data belonging to the term/pathway.
  }
}

\references{

See \link[fgsea]{fgsea} for full details of GSEA implentation in R. \cr

See \url{https://www.gsea-msigdb.org/gsea/index.jsp} for GSEA documentation. \cr

If you use the \code{gseaAnalysis} function in your analysis, please cite: \cr

A. Subramanian, P. Tamayo, V.K. Mootha, S. Mukherjee, B.L. Ebert, M.A. Gillette, A. Paulovich, S.L. Pomeroy, T.R. Golub, E.S. Lander, & J.P. Mesirov. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles, Proc. Natl. Acad. Sci. U.S.A. 102 (43) 15545-15550, (2005). \cr

Mootha, V., Lindgren, C., Eriksson, KF. et al. PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nat Genet 34, 267–273 (2003). 
}

\seealso{
\link[GSEABase]{GSEABase} \cr
\link[msigdb]{msigdb} \cr
\link[fgsea]{fgsea} \cr
\link[postNet]{postNetStart} \cr
\link[postNet]{slopeFilt} \cr
\link[postNet]{ptn_GSEA} \cr
\link[postNet]{gseaPlot} \cr
}

\examples{

# ------------------------------------------------------------------
# Example 1: Running GSEA with a custom gene set
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

# Create example custom gene sets for GSEA:
set1 <- sample(myGenes[[1]], 100)
set2 <- sample(myGenes[[2]], 100)
inSet <- list(Set1 = set1, Set2 = set2)

# Run GSEA on custom gene sets:
 ptn <- gseaAnalysis(ptn = ptn,
                    geneSet = inSet,
                    name = 'example')

# Extract the signficant enrichment results from the postNetData object:
 gseaOut <- ptn_GSEA(ptn, threshold=0.05)
 str(gseaOut)
 
# Plot GSEA results:
 gseaPlot(ptn = ptn,
          termNames = gseaOut$Term[1])
          
\dontrun{

# ------------------------------------------------------------------
# Example 2: Running GSEA with MSigDB collections
# ------------------------------------------------------------------

# If using a GSEA collection from MSigDB, check the available versions:
 version <- msigdb::getMsigdbVersions()
 str(version)

# Retrieve the MSigDB data for "human":
 msigdbOut <- msigdb::getMsigdb(org = 'hs', id = 'SYM', version = version)

# Check the available collections or subcollections:
 msigdb::listCollections(msigdbOut)
 msigdb::listSubCollections(msigdbOut)
 
# Initialize a postNetData object:
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

# Run GSEA on the C5 collection with GO:BP and specific terms:
 ptn <- gseaAnalysis(ptn = ptn,
                    collection = 'c5',
                    subcollection = 'GO:BP',
                    subsetNames = c('GOBP_CELL-CELL_SIGNALING_BY_WNT',
                    'GOBP_ENDOCYTOSIS'),
                    name = 'c5_example')

# ------------------------------------------------------------------
# Example 3: Running GSEA requiring anota2seq analysis slope filtering
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