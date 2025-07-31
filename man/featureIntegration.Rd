\name{featureIntegration}
\alias{featureIntegration}
\title{Feature integration with post-transcriptional network modelling
}

\description{
}

\usage{
featureIntegration(ptn, 
                  features, 
                  lmfeatGroup = NULL, 
                  lmfeatGroupColour = NULL, 
                  analysis_type, 
                  regOnly = TRUE, 
                  allFeat = FALSE, 
                  useCorel = TRUE, 
                  covarFilt = 20, 
                  NetModelSel = "omnibus", 
                  comparisons = NULL, 
                  fdrUni = 0.05, 
                  stepP = 0.05, 
                  pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{features}{A named list of numeric vectors corresponding to features that will be included in modelling of regulatory effects. Each vector element in the list must have names corresponding to gene IDs. These vectors quantifying features can be the outputs of previous steps of the analysis, including from \code{contentAnalysis}, \code{lengthAnalysis}, etc.. Custom inputs can also be included if they are named numeric vectors with geneIDs matching those in the \code{postNetData} object. Finally, lists of geneIDs without numeric values, for example genes regulated downstream of a pathway, can be converted to binary feature vectors using the \link[postNet]{signCalc} function.  
}
  \item{lmfeatGroup}{If \code{analysis_type = "lm"}, optionally provide an character vector specifying how \code{features} will be grouped in the network plot. Groups of features will be circled, with circle colours corresponding to the groups specified. For example, if the input list for \code{features} includes 5'UTR length, uORFs, enriched codons, 3'UTR length, and a 3'UTR motif, a grouping vector could be: \code{c("UTR5", "UTR5", "CDS", "UTR3", "UTR3")}, to organize the layout of features in the network plot according to regions of the mRNA sequence. Grouping vectors must be the same length as the input \code{features} list. Default is \code{NULL}. 
}
  \item{lmfeatGroupColour}{If \code{analysis_type = "lm"}, optionally provide a named character vector specifying colours for the different groups defined with \code{lmfeatGroup} that will be used in the network plot. Colour vectors must be the same length as the number of groups in \code{lmfeatGroup}, and be named according to group. For example, \code{c("red", "blue", "green")} with names \code{c("UTR5", "CDS", "UTR3")}. Default is \code{NULL}.
}
  \item{analysis_type}{A string specifying the method for the analysis. The options are \code{"lm"} to perform stepwise linear regression modelling, or \code{"rf"} to perform a random forest-based analysis. 
}
  \item{regOnly}{Logical specifying whether modelling should be performed using only regulated genes (i.e., those included in \code{geneList}), or with all genes in the datase. If \code{regOnly = TRUE}, only regulated genes will be used, if \code{regOnly = FALSE}, all genes will be considered. Default is \code{TRUE}. 
}
  \item{allFeat}{If \code{analysis_type = "lm"}, provide a logical specifying which input features should be including in the stepwise linear regression modelling. If \code{allFeat = TRUE}, all input features will be included, if \code{FALSE}, only features passing the FDR threshold specified with \code{fdrUni} in univariate models will be included. Default is \code{FALSE}.
}
  \item{useCorel}{
}
  \item{covarFilt}{Minimum threshold for covariance between features to be dispayed for the network plot. The default is 20.
}
  \item{NetModelSel}{Values for the variance explained for the network plot. The default is "Omnibus". Aternatively "Adjusted" or "Univariate" can be provided. 
}
  \item{comparisons}{If \code{regOnly = TRUE}, provide a list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}. If both  \code{regulation} and \code{geneList} are provided, \code{geneList} is appended to the end.
}
  \item{fdrUni}{
}
  \item{stepP}{
}
  \item{pdfName}{Name of the output pdf files.
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
