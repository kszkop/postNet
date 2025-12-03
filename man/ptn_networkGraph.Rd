\docType{methods}
\name{ptn_networkGraph}
\alias{ptn_networkGraph,postNetData-method}
\title{Accessor for the \code{networkGraph} slot of a \code{postNetData} object
}

\description{Retrieve the network (an S3 \link[igraph]{igraph} object) resulting from stepwise regression modelling from the \code{postNetData} object.
}

\usage{
\S4method{ptn_networkGraph}{postNetData}(x, comparison)
}

\arguments{
  \item{x}{A \linkS4class{postNetData} object.
}
  \item{comparison}{An integer specifying which comparison to select results from. For example, \code{1} would
  select the first comparison specified in the \code{comparisons} parameter of the \link[postNet]{featureIntegration} function.
}
}

\value{An S3 \link[igraph]{igraph} object based on the results of stepwise regression modelling with \link[postNet]{featureIntegration}. 
}

\seealso{
\link[postNet]{featureIntegration} \cr
\link[igraph]{igraph}
}

\examples{
# load and create example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn
 
# Prepare the list of pre-calculated features to be used in feature integration modelling 
 myFeatures <- postNetExample$features
 str(features)

# Run feature integration modelling using stepwise regression
 ptn <- featureIntegration(ptn = ptn, 
                          features = myFeatures, 
                          pdfName = "omnibus", 
                          regOnly = T, 
                          allFeat = F,
                          analysis_type = "lm",
                          covarFilt = 10,
                          comparisons = list(c(1,2)),
                          NetModelSel = "omnibus")
                          
# Retrieve the network graph of the results of the stepwise regression omnibus model
 networkGraph <- ptn_networkGraph(ptn, comparison = 1) 
 str(networkGraph)
}
