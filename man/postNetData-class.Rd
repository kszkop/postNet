\name{postNetData-class}
\Rdversion{1.1}
\docType{class}
\alias{postNetData-class}
\alias{ptn_codonSelection,postNetData-method}
\alias{ptn_features,postNetData-method}
\alias{ptn_motifgeneList,postNetData-method}
\alias{ptn_motifSelection,postNetData-method}
\alias{ptn_networkGraph,postNetData-method}

\title{Class \code{"postNetData"}}
\description{
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("postNetData", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{species}:}{Object of class \code{"characterOrNULL"} ~~ }
    \item{\code{version}:}{Object of class \code{"characterOrNULL"} ~~ }
    \item{\code{selection}:}{Object of class \code{"character"} ~~ }
    \item{\code{annot}:}{Object of class \code{"postNetAnnot"} ~~ }
    \item{\code{dataIn}:}{Object of class \code{"postNetDataIn"} ~~ }
    \item{\code{features}:}{Object of class \code{"dataframeOrNULL"} ~~ }
    \item{\code{analysis}:}{Object of class \code{"postNetAnalysis"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{ptn_background}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_codonAnalysis}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_codonSelection}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_colours}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_dataIn}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_effect}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_features}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_geneID}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_geneList}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_id}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_motifgeneList}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_motifSelection}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_networkGraph}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_selection}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_sequences}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_species}{\code{signature(x = "postNetData")}: ... }
    \item{ptn_version}{\code{signature(x = "postNetData")}: ... }
	 }
}
\references{
}
\author{
}
\note{
}

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("postNetData")
}
\keyword{classes}
