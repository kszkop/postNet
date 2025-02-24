\name{gseaAnalysis}
\alias{gseaAnalysis}
\title{Perform Gene Set Enrichment Analysis (GSEA) with a \code{postNetData} object 
}

\description{
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
  \item{genesSlopeFiltOut}{If using an \pkg{anota2seq} analysis with gene sets from the "translation" or "buffering" 
  regulatory modes, optionally provide a character vector with gene IDs to be excluded from GSEA that is the
  output of the \code{\link[postNet]{slopeFilt}} function. Although optional, this step is strongly recommended. This 
  filtering is not necessary if using "mRNAAbundance", "totalmRNA", or "translatedmRNA". Default is \code{NULL}.  
}
  \item{collection}{A character selecting the MSigDB collection of interest, for example:  \code{'c8'}, or \code{'h'}.
  Although possible, it is not recommended to run all collections in the same analysis. Available options for
  collections can be found at \url{https://www.gsea-msigdb.org/gsea/msigdb/index.jsp}, or by using the \pkg{msigdb}
  package (see example below). Currently, collections are only avaible for human and mouse. 
}
  \item{subcollection}{Optionally, provide a character selecting the MSigDB subcollection of interest. 
  For example: \code{"GO:BP"}, or \code{"CP:REACTOME"}. Available options for subcollections can be found at 
  \url{https://www.gsea-msigdb.org/gsea/msigdb/index.jsp}, or by using the \pkg{msigdb} package (see example below).
  Currently, collections are only avaible for human and mouse. Default is \code{NULL}. 
}
  \item{subsetNames}{
}
  \item{geneSet}{
}
  \item{maxSize}{
}
  \item{minSize}{
}
  \item{name}{
}
}

\details{
}

\value{A data.frame of GSEA results that is stored in the GSEA slot of the \code{postNetData} object, and written to a tab delimited text file. Each row of the data.frame corresponds to a term, with the following columns (see \code{\link[fgsea]{fgseaMultilevel}} from the \pkg{fgsea} package:

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

See \code{\link[fgsea]{fgsea}} for full details of GSEA implentation in R.

See \url{https://www.gsea-msigdb.org/gsea/index.jsp} for GSEA documentation. \cr

A. Subramanian, P. Tamayo, V.K. Mootha, S. Mukherjee, B.L. Ebert, M.A. Gillette, A. Paulovich, S.L. Pomeroy, T.R. Golub, E.S. Lander, & J.P. Mesirov. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles, Proc. Natl. Acad. Sci. U.S.A. 102 (43) 15545-15550, (2005). \cr

Mootha, V., Lindgren, C., Eriksson, KF. et al. PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nat Genet 34, 267–273 (2003). 
}

\note{
}


\seealso{
\code{\link[GSEABase]{GSEABase}} \cr
\code{\link[msigdb]{msigdb}} \cr
\code{\link[fgsea]{fgsea}} \cr
\code{\link[postNet]{slopeFilt}} \cr
\code{\link[postNet]{ptn_GSEA}} \cr
\code{\link[postNet]{gseaPlot}} \cr
}

\examples{
}
