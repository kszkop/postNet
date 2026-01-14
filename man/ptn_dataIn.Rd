\docType{methods}
\name{ptn_dataIn}
\alias{ptn_dataIn,postNetData-method}
\title{Accessor for the "background", "geneList", "effect", and "colours" slots of a \code{postNetData} object
}

\description{Retrieve the "background", "geneList", "effect", and "colours" slots from a \code{postNetData} object. If \code{ads} has been provided, these are retrieved from the \code{Anota2seqDataSet} based on the \code{regulation}, \code{contrast}, \code{regulationGen}, and \code{contrastSel} arguments of the \code{\link[postNet]{postNetStart}} function. If not using \code{ads}, these inputs are compiled from the \code{geneList}, \code{geneListcolours}, \code{effectMeasure}, and optionally the \code{customBg} arguments.
}

\usage{
\S4method{ptn_dataIn}{postNetData}(ptn)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }
}

\value{Returns an S4 object of class \code{postNetDataIn} with the following slots:
  
  \itemize{ 
  \item \strong{background}: a character vector of gene/transcript IDs corresponding to the background used for
  statistical comparisons against the gene sets of interest in \code{geneList}. This background is either taken
  automatically from the \code{Anota2seqDataSet} if using \code{ads}, or if not using \code{ads} it is created from 
  a user-provided \code{geneList} or provided directly with \code{customBg}.
  
  \item \strong{geneList}: a list of one or more character vectors corresponding to the gene sets of interest to be
  used in the \pkg{postNet} analysis. If using \code{ads}, gene sets corresponding to \pkg{anota2seq} regulatory modes
  are selected using the \code{regulation} argument. If not using \code{ads}, the list of gene sets is provided using
  \code{geneList}.
  
  \item \strong{effect}: a numeric vector corresponding to the regulation effect measurement. The elements of the vector
  are named with the gene/transcript IDs present in either \code{ads}, \code{geneList}, or \code{customBg}. If using
  \code{ads}, the regulatory effect measurement is selected by \code{regulationGen}. If using \code{geneList}, it is
  provided with \code{effectMeasure}. 
  
  \item \strong{colours}: a character vector of the same length as \code{geneList} specifying the colours to be
  used for plotting. If using \code{geneList}, colours are user-provided. If using ads, the default colours for each 
  regulatory mode from \pkg{anota2seq} will be used automatically.

  }
}
  
\examples{

# load and create example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Access the input data (background, geneList, effect, colours) stored in the postNetData object:
  myDataIn <- ptn_dataIn(ptn)
  str(myDataIn)
  
}

\seealso{
\code{\link[postNet]{postNetStart}} for more information on required inputs for postNet analysis. \cr
\code{\link[postNet]{postNetData-class}}
}
