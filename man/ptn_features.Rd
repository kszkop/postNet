\docType{methods}
\name{ptn_features}
\alias{ptn_features,postNetData-method}
\title{Accessor for the \code{features} slot of a \code{postNetData} object
}

\description{Retrieve the input features used in stepwise regression and random forest modelling from the \code{postNetData} object.
}

\usage{
\S4method{ptn_features}{postNetData}(x)
}

\arguments{
  \item{x}{A \linkS4class{postNetData} object.
}
}

\value{Returns a \code{data.frame} of features and/or gene signatures (see \link[postNet]{signCalc}) quantified during a postNet analysis and used as input for \link[postNet]{featureIntegration}. Rows correspond to genes, and columns correspond to features.  
}

\seealso{
\link[postNet]{featureIntegration}
\link[postNet]{signCalc}
}

\examples{
}
