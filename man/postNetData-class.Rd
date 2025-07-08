\name{postNetData-class}
\docType{class}
\alias{postNetData,postNetData-class}

\title{S4 Class \code{postNetData}: Core data container for postNet analyses
}

\description{The \code{postNetData} class stores the core data structure used throughout the postNet analysis. It contains input gene lists, sequence annotations and source information, and results from codon, motif, GSEA, GAGE, Go term, and feature integration analyses. This object is created by the \code{\link[postNet]{postNetStart}} function and serves as the input/output container for most downstream functions.
}

\section{Objects from the Class}{
Objects of this class are typically created using the constructor function \code{\link[postNet]{postNetStart}}. Alternatively, they can be created manually using \code{new("postNetData", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{species}:}{Object of class \code{"characterOrNULL"} Specifies the species associated with the sequence annotations and input gene lists, e.g. \code{"human"}, \code{"mouse"}. Used to load species-specific reference sequences.
 }
    \item{\code{version}:}{Object of class \code{"characterOrNULL"} Version identifier for the RefSeq sequence annotations loaded.
}
    \item{\code{selection}:}{Object of class \code{"character"} Indicates the method used to select representative transcript isoforms, e.g. \code{"random"}, \code{"shortest"}, or \code{"longest"}.
}
    \item{\code{seed}:}{Object of class \code{"ANY"} Stores the integer used to set the seed for random sampling if \code{selection = "random"}
}
    \item{\code{annot}:}{Object of class \code{"postNetAnnot"} Stores the annotated sequences for each mRNA transcript region, \code{"UTR5"}, \code{"CDS"}, \code{"CCDS"}, and \code{"UTR3"}.
}
    \item{\code{dataIn}:}{Object of class \code{"postNetDataIn"} Contains the user-defined input gene lists of interest, background gene set, regulatory effect measurement, and selected colours for plotting.
}
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
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("postNetData")
}
\keyword{classes}
