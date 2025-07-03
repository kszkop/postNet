\name{featureIntegration}
\alias{featureIntegration}
\title{Post-transcriptional network modelling and feature integration
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
                  pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{features}{Named list of numeric vectors corresponding to features that will be included in modelling the explain variance in the regulatory effect measurement. Each vector element in the list must have names corresponding to gene IDs. These vectors quantifying features can be the outputs of previous steps of the analysis, including from \code{contentAnalysis}, \code{lengthAnalysis}, etc.. Custom inputs can also be included if they are named numeric vectors with geneIDs matching those in the \code{postNetData} object. Finally, lists of geneIDs without numeric values, for example genes regulated downstream of a pathway, can be converted to binary feature vectors using the \link[postNet]{signCalc} function.      
}
  \item{lmfeatGroup}{An optional character vector specifying how \item{features} will be grouped in the network plot. Groups of features will be circled, with circle colours corresponding to the groups specified. For example, if the input list for \item{features} includes 5'UTR length, uORFs, enriched codons, 3'UTR length, and a 3'UTR motif, a grouping vector could be: \code{c("UTR5", "UTR5", "CDS", "UTR3", "UTR3")}, to organize the layout of features in the network plot according to regions of the mRNA sequence. Grouping vectors must be the same length as the input \code{features} list. Default is \code{NULL}.      
}
  \item{lmfeatGroupColour}{
}
  \item{analysis_type}{Method for the analysis. Choose "lm" for stepwise linear model or "rf" random forest based analysis
}
  \item{regOnly}{Logical, whether only regulated genes should be considered for the analysis or all expressed. The default is TRUE. 
}
  \item{allFeat}{If the analysis_type is "lm". Please select whether all features from input should be selected for the stepwise linear model or only these 
  that passed 0.05 fdr threshold in the step of univariate models. The default is TRUE.
}
  \item{useCorel}{
}
  \item{covarFilt}{Minimum threshold for covariance between features to be dispayed for the network plot. The default is 20. 
}
  \item{NetModelSel}{Values for the variance explained for the network plot. The default is "Omnibus". Aternatively "Adjusted" or "Univariate" can be provided. 
}
  \item{comparisons}{Relevant for linear model only (analysis_type="lm") and if regulated genes are only considered (regOnly=TRUE). Indicate which regulation/geneList pairs should be compared. It should be provided as a list of vectors with pairs of indexes in regulation/geneList (for example list(c(0,2),c(2,4)) would indicate that the background (0) would be compared to the second and the second to the fourth mode of regulation/geneList. If both regulation and geneList are provided, geneList is appended to the end.
  
  A list of numeric vectors indicating which gene sets defined in the \code{regulation} or \code{geneList} slots
  of the \code{postNetData} object should be used in pairwise comparisons. The numeric vectors should correspond to the indexes of
  the gene sets, where 0 denotes the background. For example, \code{list(c(0,2),c(1,2))} indicates that the second gene set in
  \code{regulation} or \code{geneList} will be compared to the background (0), and the first and second gene sets will be compared
  to each other. The default is \code{NULL}.
  
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
