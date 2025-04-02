\name{goAnalysis}
\alias{goAnalysis}
\title{Perform Gene Ontology Analysis on gene lists of interest with a \code{postNetData} object
}

\description{The \code{goAnalysis} function implements \link[clusterProfiler]{clusterProfiler} to perform Gene Ontology 
enrichment analysis using the regulatory effect measurement contained in a \code{postNetData} object. These gene sets 
are specified by the \code{regulationGen}, or \code{effectMeasure} parameters in the \link[postNet]{postNetStart} function. 

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
  Terms with fewer genes than the threshold specified will be discarded, regardless of significance. Default is \code{10}.
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

\details{
}
\value{
}
\references{
If you use the \code{goAnalysis} function in your analysis, please cite: \cr

S Xu, E Hu, Y Cai, Z Xie, X Luo, L Zhan, W Tang, Q Wang, B Liu, R Wang, W Xie, T Wu, L Xie, G Yu.
  Using clusterProfiler to characterize multiomics data. Nature Protocols. 2024, 19(11):3292-3320
}

\note{
}

\seealso{
\link[clusterProfiler]{clusterProfiler}
\link[postNet]{postNetStart} \cr
\link[postNet]{slopeFilt} \cr
\link[postNet]{ptn_GO} \cr
\link[postNet]{goDotPlot} \cr

}
\examples{
}
