\name{ptn_GO}
\alias{ptn_GO}
\title{Accessor for the "GO" slot of a \code{postNetData} object
}

\description{Retrieve the \code{"GO"} slot from a \code{postNetData} object, which holds the results of GO term enrichment analysis created by 
the \link[postNet]{goAnalysis} function.
}

\usage{
ptn_GO(ptn, 
      category, 
      geneList, 
      threshold)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{category}{A character vector specifying one or more Gene Ontology categories to extract results for. 
  The options are \code{'BP'}, \code{'CC'}, \code{'MF'}, and \code{'KEGG'}. For example, \code{category = c('BP','MF')}.
}
  \item{geneList}{A string specifying the name of the gene set of interest to extract results for. This must match one of the names 
  of the elements of \code{geneList} in the \code{postNetData} object. 
}
  \item{threshold}{A numeric value providing the FDR threshold to be used to select results from the GO term analysis.
}
}

\value{A data.frame of GO term enrichment analysis results that is stored in the GO slot of the \code{postNetData} object, and written to an Excel file. Each row of the data.frame corresponds to a term, with the following columns:

 \itemize{ 
  \item \strong{ID}: The ID of the term or pathway tested.
  
  \item \strong{Description}: The name of the term. Note that if the \code{category} used in the analysis was "KEGG", there are two additional
  columns with descriptive information called "category" and "subcategory".
  
  \item \strong{Count}: The number of genes in the set belonging to the term or pathway
  
  \item \strong{Size}: The size of the pathway after removing genes not present in the input.  
  
  \item \strong{pvalue}: The enrichment p-value.  
  
  \item \strong{adjusted_pvalue}: The BH-adjusted p-value.
  
  \item \strong{geneID}: The gene IDs in the input data belonging to the term/pathway.
  }
}

\references{
See \link[clusterProfiler]{clusterProfiler} for full details of enrichment analysis implementation in R. \cr

If you use the \code{goAnalysis} function in your analysis, please cite: \cr

S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang, W Xie, T Wu, L Xie, G Yu.
  Using clusterProfiler to characterize multiomics data. Nature Protocols. 2024, 19(11):3292-3320
}

\seealso{
\link[clusterProfiler]{clusterProfiler} \cr
\link[postNet]{goAnalysis} \cr
\link[postNet]{postNetStart} \cr
\link[postNet]{slopeFilt} \cr
\link[postNet]{goDotplot} \cr
}

\examples{
# load example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Run GO term analysis:
 ptn <- goAnalysis(ptn = ptn, 
                  category=c("BP"),
                  name = "example")

# Extract the significant enrichment results from the postNetData object:
 goOut <- ptn_GO(ptn, 
                category = "BP", 
                geneList = "translationUp", 
                threshold = 0.05)
 str(goOut)
}
