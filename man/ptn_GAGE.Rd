\docType{function}
\name{ptn_GAGE}
\alias{ptn_GAGE}
\title{Accessor for the "GAGE" slot of a \code{postNetData} object
}

\description{Retrieve the "GAGE" slot from a \code{postNetData} object, which holds the results of GAGE analysis created by the \link[postNet]{gageAnalysis} function.
}

\usage{
ptn_GAGE(ptn,
        category, 
        direction, 
        threshold)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }
  \item{category}{A character vector specifying one or more gene ontology categories to be included in the analysis. 
  The options are \code{"BP"}, \code{"CC"}, \code{"MF"}, and \code{"KEGG"}. For example, \code{category = c("BP", "MF")}.
  }
  \item{direction}{A character indicating the directionality of regulation of the gene sets. Can be either
  \code{"greater"}, or \code{"less"}.
  }
  \item{threshold}{A numeric value providing the FDR threshold to be used to select results from the GAGE analysis. 
  }
}

\value{The value returned by the \code{ptn_GAGE} function is a data.frame with the following columns (see the \link[gage]{gage} package for full details):

 \itemize{ 
  \item \strong{p.geomean}: The geometric mean of the individual p-values from multiple single array-based gene set 
  tests.
  
  \item \strong{stat.mean}: The mean of the individual statistics from multiple single array-based gene set tests. The
  value denotes the magnitude of the gene-set level changes, and the sign denotes the direction of the changes.
  
  \item \strong{p.val}: The global p-value or summary of the individual p-values from multiple single array-based gene
  set tests.
  
  \item \strong{q.val}: The BH adjusted p-value, as implemented by the \pkg{multtest} package.
  
  \item \strong{set.size}: The number of genes included in the gene set.  
  
  \item \strong{Genes}: The gene IDs in the input data belonging to the gene ontology gene set.
  }
}

\references{
See the \link[gage]{gage} package for details. \cr

Luo, W., Friedman, M., Shedden K., Hankenson, K. and Woolf, P GAGE: Generally Applicable Gene Set Enrichment for Pathways Analysis. BMC Bioinformatics 2009, 10:161.

}

\seealso{
\link[gage]{gage} \cr
\link[postNet]{gageAnalysis} \cr
\link[postNet]{postNetStart} \cr
}

\examples{
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
                    threshold = 1)
 str(gageOut)
}