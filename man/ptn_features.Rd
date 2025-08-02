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

\value{Returns a \code{data.frame} of features (or gene signatures, see \link[postNet]{signCalc}) quantified during a postNet analysis and used as input for \link[postNet]{featureIntegration}. Rows correspond to genes, and columns correspond to features.  
}

\seealso{
\link[postNet]{featureIntegration}
\link[postNet]{signCalc}
}

\examples{

# load and create example data:
 data(postNetExample)

 myGenes <- postNetExample$geneList
 myBg <- postNetExample$background
 myEffect <- postNetExample$effect

# Initialize the postNetData object with custom gene lists:
 ptn <- postNetStart(
   geneList = myGenes,
   geneListcolours = c("#7FB7BE","#DB7F67"),
   customBg = myBg,
   effectMeasure = myEffect,
   selection = "random",
   source = "load", 
   species = "human"
 )
 
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
                          covarFilt = 20,
                          comparisons = list(c(1,2)),
                          NetModelSel = "omnibus")
                          
# Retrieve a data.frame summary of the features used in modelling
 featureInput <- ptn_features(ptn)
 str(featureInput)
}
