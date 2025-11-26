\docType{methods}
\name{ptn_background}
\alias{ptn_background,postNetData-method}
\title{Accessor for the gene set used as background from the "background" slot of a \code{postNetData} object
}
\description{The "background" slot holds a character vector of gene IDs used as the background in statistical comparisons performed in \pkg{postNet} analyses.
}
\usage{
\S4method{ptn_background}{postNetData}(x)
}

\arguments{
  \item{x}{a \code{postNetData} object.
  }
}

\value{A character vector containing the gene IDs for the set of genes used as background.
}

\seealso{
\code{\link[postNet]{postNetStart}} for details on background gene sets.
}

\examples{

# load and create example data:
  postNetExample <- postNetParcel::postNetExample
  ptn <- postNetExample$ptn

# Access the set of background genes used in postNet analyses
  myBg <- ptn_background(ptn)
  str(myBg)

}
