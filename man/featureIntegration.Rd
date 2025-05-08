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
  \item{ptn}{
}
  \item{features}{Named list of features. Each entry in the list is a vector with "quantified" feature at each step of the analysis.
}
  \item{lmfeatGroup}{
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
  \item{comparisons}{Relevant for linear model only (analysis_type="lm") and if regulated genes are only considered (regOnly=TRUE). Indicate which 
  regulation/geneList pairs should be compared. It should be provided as a list of vectors with pairs of indexes in regulation/geneList (for example list(c(0,2),c(2,4)) 
  would indicate that the background (0) would be compared to the second and the second to the fourth mode of regulation/geneList. If both regulation and geneList are 
  provided, geneList is appended to the end.
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
