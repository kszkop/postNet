\name{goDotplot}
\alias{goDotplot}
\title{Produce dot plots with results of GO term enrichment analysis
}

\description{The \code{goDotplot} function can be used to produce dot plots visualizing the results of GO term enrichment analyses
performed using \link[postNet]{goAnalysis}. The input is a \code{postNetData} object, with several options for plotting.
}

\usage{
goDotplot(ptn, 
          category, 
          pool = TRUE, 
          termSel = NULL, 
          nCategories = 10, 
          size = "Count", 
          pdfName = NULL)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
}
  \item{category}{A character vector specifying one or more Gene Ontology categories to plot results for. 
  The options are \code{'BP'}, \code{'CC'}, \code{'MF'}, and \code{'KEGG'}. For example, \code{category = c('BP','MF')}
}
  \item{pool}{Logical indicating if results for all gene sets of interest should be pooled into a single plot. For example,
  if \code{code}, results for "translationUp" and "translationDown" will be displayed together, automatically using the colours 
  stored in the \code{postNetData} object. If \code{FALSE}, a separate plot will be generated for each gene set of interest. 
  Default is \code{TRUE}.
}
  \item{termSel}{A character vector can be supplied to specify specific terms to be plotted. For example, \code{c('GO:0008380','GO:0008380')},
  or \code{c('hsa03013')} is using KEGG terms. Default is \code{NULL}.  
}
  \item{nCategories}{An integer specifying the maximum number of terms to include on each plot. Default is \code{10}. 
}
  \item{size}{A string specifying the metric used to determine the size of the dots. The options are \code{"Count"}, which is the
  number of genes from the gene set of interest included in each specific term, or \code{"geneRatio"} which is the "Count" divided
  by the total number of genes in the term. Default is \code{"Count"}. 
}
  \item{pdfName}{Name to be appended to output PDF files. Default is \code{NULL}.
}
}

\details{
}

\value{
}

\references{
}

\seealso{
}

\examples{
}
