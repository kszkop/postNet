\name{plotFeaturesMap}
\alias{plotFeaturesMap}
\title{Explore mRNA features and regulatory effects using UMAP visualizations 
}

\description{
}

\usage{
plotFeaturesMap(ptn, 
                regOnly = TRUE, 
                comparisons = NULL, 
                featSel, 
                remBinary = TRUE, 
                featCol = NULL, 
                scaled = FALSE, 
                centered = TRUE, 
                remExtreme = NULL, 
                pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{regOnly}{Logical specifying whether UMAPs should include only regulated genes (i.e., those included in \code{geneList}), or with all genes in the dataset. If \code{regOnly = TRUE}, only regulated genes will be used, if \code{regOnly = FALSE}, all genes will be considered. Default is \code{TRUE}. 
}
  \item{comparisons}{If \code{regOnly = TRUE}, provide a list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}. If both  \code{regulation} and \code{geneList} are provided, \code{geneList} is appended to the end so numbering will continue sequentially.
}
  \item{featSel}{
}
  \item{remBinary}{Default is \code{TRUE}.
}
  \item{featCol}{Default is \code{NULL}.
}
  \item{scaled}{Default is \code{FALSE}.
}
  \item{centered}{Default is \code{TRUE}.
}
  \item{remExtreme}{Default is \code{NULL}.
}
  \item{pdfName}{Name to be appended to output PDF files.
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
