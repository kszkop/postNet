\name{ptn_GSEA}
\alias{ptn_GSEA}
\title{Accessor for the "GSEA" slot of a \code{postNetData} object
}

\description{Retrieve the "GSEA" slot from a \code{postNetData} object, which holds the results of GSEA created by the \link[postNet]{gseaAnalysis} function.
}

\usage{
ptn_GSEA(ptn, threshold)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
 }
  \item{threshold}{A numeric value providing the FDR threshold to be used to select results from the GSEA. 
 }
}

\value{A data.frame of GSEA results that is stored in the GSEA slot of the \code{postNetData} object, and written to a tab-delimited text file. Each row of the data.frame corresponds to a term, with the following columns (also see \code{\link[fgsea]{fgseaMultilevel}} from the \pkg{fgsea} package):

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
See \code{\link[fgsea]{fgsea}} for full details of GSEA implementation in R. \cr

See \url{https://www.gsea-msigdb.org/gsea/index.jsp} for GSEA documentation. \cr

A. Subramanian, P. Tamayo, V.K. Mootha, S. Mukherjee, B.L. Ebert, M.A. Gillette, A. Paulovich, S.L. Pomeroy, T.R. Golub, E.S. Lander, & J.P. Mesirov. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles, Proc. Natl. Acad. Sci. U.S.A. 102 (43) 15545-15550, (2005). \cr

Mootha, V., Lindgren, C., Eriksson, KF. et al. PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nat Genet 34, 267–273 (2003).
}

\seealso{
\code{\link[fgsea]{fgsea}} \cr
\code{\link[postNet]{gseaAnalysis}} \cr
\code{\link[postNet]{gseaPlot}}
}

\examples{
# load example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Create example custom gene sets for GSEA:
 inSet <- list(set=sample(unlist(postNetExample$geneList[[1]]), 10))

# Run GSEA on custom gene sets:
 ptn <- gseaAnalysis(ptn = ptn,
                    geneSet = inSet,
                    name = 'example')

# Extract enrichment results from the postNetData object with an FDR < 0.05:
 gseaOut <- ptn_GSEA(ptn, threshold=0.05)
}

