\name{gseaPlot}
\alias{gseaPlot}
\title{Plot enriched terms from GSEA
}

\description{The \code{gseaPlot} function generates visualizations of the results of GSEA. The rank is plotted on the x-axis against the enrichment score on the y-axis. The positions of genes belonging to the term or pathway of interest are indicated in the ranking by red bars.   
}

\usage{
gseaPlot(ptn, 
        termNames, 
        genesSlopeFiltOut = NULL, 
        gseaParam = 1, 
        ticksSize = 0.3, 
        pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }
  \item{termNames}{Character vector specifying the names of the terms to be plotted.
  }
  \item{genesSlopeFiltOut}{If using an \code{Anota2seqDataSet} with analysis results for the "translation" or
  "buffering" regulatory modes, optionally provide a character vector with gene IDs to be excluded from GSEA that 
  is the output of the \code{\link[postNet]{slopeFilt}} function. Although optional, this step is strongly 
  recommended. This filtering is not necessary if using "mRNAAbundance", "totalmRNA", or "translatedmRNA". 
  Default is \code{NULL}.
  }
  \item{gseaParam}{The GSEA parameter used to adjust the displayed statistic values. Can be either 0 or 1. 
  Default is \code{1}. See \code{\link[fgsea]{fgseaMultilevel}} from the \pkg{fgsea} package for details.
  }
  \item{ticksSize}{Number to adjust the thickness of the lines indicating the position of genes in the ranking. Default 
  is \code{0.3}. 
  }
  \item{pdfName}{Optionally, provide a name for the output PDF file.
  }
}

\value{No value is returned. A graphical output is generated as a PDF file.
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
\code{\link[postNet]{ptn_GSEA}}
}

\examples{
# load example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Create example custom gene sets for GSEA:
 inSet <- list(set = sample(unlist(postNetExample$geneList[[1]]), 10))

# Run GSEA on custom gene sets:
 ptn <- gseaAnalysis(ptn = ptn,
                    geneSet = inSet,
                    name = "example")

# Extract the significant enrichment results from the postNetData object:
 gseaOut <- ptn_GSEA(ptn, threshold = 0.05)

# Plot GSEA results:
 gseaPlot(ptn = ptn,
          termNames = gseaOut$Term[1],
          pdfName = "my_GSEAplot")

}

