\docType{methods}
\name{ptn_features}
\alias{ptn_features,postNetData-method}
\title{Accessor for the \code{features} slot of a \code{postNetData} object
}

\description{Retrieve the input features used in stepwise regression and random forest modelling from the \code{postNetData} object.
}

\usage{
\S4method{ptn_features}{postNetData}(ptn)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
}

\value{Returns a \code{data.frame} of features (or gene signatures, see \link[postNet]{signCalc}) quantified during a postNet analysis and used as input for \link[postNet]{featureIntegration}. Rows correspond to genes, and columns correspond to features.  
}

\seealso{
\link[postNet]{featureIntegration} \cr
\link[postNet]{signCalc}
}

\examples{

# load and create example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn
 
# Prepare the list of pre-calculated features to be used in feature integration modelling:
 myFeatures <- postNetExample$features
 str(myFeatures)

# Run feature integration modelling using forward stepwise regression:
 ptn <- featureIntegration(ptn = ptn, 
                          features = myFeatures, 
                          pdfName = "omnibus", 
                          regOnly = TRUE, 
                          allFeat = FALSE,
                          analysis_type = "lm",
                          covarFilt = 20,
                          comparisons = list(c(1,2)),
                          NetModelSel = "omnibus")
                          
# Retrieve a data.frame summarizing the features for each gene used in modelling:
 featureInput <- ptn_features(ptn)
 str(featureInput)
}
