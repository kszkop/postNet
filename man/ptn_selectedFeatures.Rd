\docType{methods}
\name{ptn_selectedFeatures}
\alias{ptn_selectedFeatures,postNetData-method}
\title{Accessor for the \code{selectedFeatures} slot of a \code{postNetData} object
}
\description{Retrieve the features explaining changes in the regulatory effect identified by stepwise regression or random forest models from the \code{postNetData} object.
}

\usage{
\S4method{ptn_selectedFeatures}{postNetData}(ptn, analysis_type)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{analysis_type}{A string specifying the method of analysis used for \link[postNet]{featureIntegration} to retrieve the selected features from. The options are \code{"lm"} to retrieve features selected by stepwise linear regression, or \code{"rf"} for features selected using random forest-based analyses.
}
}

\details{The relationships between selected features and changes in the regulatory effect can be further explored using the \link[postNet]{plotFeaturesMap} function, which uses UMAPs to visualize regulatory effects and features across genes. 

The selected features can also be quantified in other gene lists or datasets and used to predict regulation using random forest classification, implemented with the \link[postNet]{rfPred} function.  
}

\value{
}

\seealso{
\link[postNet]{featureIntegration} \cr
\link[postNet]{plotFeaturesMap} \cr
\link[postNet]{rfPred} 
}

\examples{
}
