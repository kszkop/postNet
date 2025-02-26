\name{goAnalysis}
\alias{goAnalysis}
\title{Perform Gene Ontology Analysis on gene lists of interest with a \code{postNetData} object
}

\description{

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
  \item{category}{A character vector specifying one or more gene ontology categories to be included in the analysis. 
  The options are \code{'BP'}, \code{'CC'}, \code{'MF'}, and \code{'KEGG'}. For example, \code{category = c('BP','MF')}.
  }
  \item{maxSize}{Integer specifying the maximal size of a gene set to test. Terms/pathways with more genes will 
  be exluded. Default is \code{500}.
  }
  \item{minSize}{Integer specifying the minimal size of a gene set to test. Terms/pathways with fewer genes will 
  be exluded. Default is \code{10}.
  }
  \item{counts}{
  }
  \item{FDR}{A numeric value providing the FDR threshold to be used to select significant enrichments. Default is
  \code{0.15}
  }
  \item{name}{
  }
}

\details{
}
\value{
}
\references{
}

\note{
}


\seealso{
}
\examples{
}
