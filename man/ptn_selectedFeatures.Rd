\docType{methods}
\name{ptn_selectedFeatures}
\alias{ptn_selectedFeatures,postNetData-method}
\title{Accessor for the \code{selectedFeatures} slot of a \code{postNetData} object
}
\description{Retrieve significant features identified by the step-wise regression model from the postNetData object
}

\usage{
\S4method{ptn_selectedFeatures}{postNetData}(ptn, analysis_type)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{analysis_type}{A string specifying the method of analysis to retrieve the selected features from. The options are \code{"lm"} to retrieve stepwise linear regression modelling results, or \code{"rf"} for features selected using random forest-based analysis.
}
}

\details{If \code{analysis_type = "lm"}
}

\value{
}

\seealso{
\link[postNet]{featureIntegration}
}

\examples{
}
