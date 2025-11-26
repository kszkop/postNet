\docType{methods}
\name{ptn_selectedFeatures}
\alias{ptn_selectedFeatures,postNetData-method}
\title{Accessor for the \code{selectedFeatures} slot of a \code{postNetData} object
}

\description{Retrieve the features explaining changes in the regulatory effect identified by stepwise regression or random forest models from the \code{postNetData} object.
}

\usage{
\S4method{ptn_selectedFeatures}{postNetData}(ptn, analysis_type, comparison)
}
\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{analysis_type}{A string specifying the method of analysis used for \link[postNet]{featureIntegration} to retrieve the selected features from. The options are \code{"lm"} to retrieve features selected by stepwise linear regression, or \code{"rf"} for features selected using random forest-based analyses.
}
  \item{comparison}{An integer specifying which comparison to select results from. For example, \code{1} would
  select the first comparison specified in the \code{comparisons} parameter of the \link[postNet]{featureIntegration} function.
}
}
\details{The relationships between selected features and changes in the regulatory effect can be further explored using the \link[postNet]{plotFeaturesMap} function, which uses UMAPs to visualize regulatory effects and features across genes. 
  
  The selected features can also be quantified in other gene lists or datasets and used to predict regulation using random forest classification, implemented with the \link[postNet]{rfPred} function.  
}

\value{If \code{analysis_type = "lm"} a named numeric vector of the features selected by stepwise regression modelling will be returned. Selected features are those passing the significance defined using the \code{stepP} parameter of the \link[postNet]{featureIntegration} function.  Numeric values correspond to the proportion of variance in the regulatory effect explained by each feature (from either omnibus, or adjusted models). 

If \code{analysis_type = "rf"}, the selected features identified by \link[Boruta]{Boruta} and used in random forest classification will be returned, with values corresponding to feature importance calculated by \link[randomForest]{randomForest}.  
}

\references{
Liaw A, Wiener M (2002). “Classification and Regression by randomForest.” R News, 2(3), 18-22. https://CRAN.R-project.org/doc/Rnews/. 

Witold R. Rudnicki M (2010). “Feature Selection with the Boruta Package.” Journal of Statistical Software, 36(11), 1–13. https://doi.org/10.18637/jss.v036.i11.
}

\seealso{
\link[postNet]{featureIntegration} \cr
\link[postNet]{plotFeaturesMap} \cr
\link[postNet]{rfPred} \cr
\link[Boruta]{Boruta} \cr
\link[randomForest]{randomForest}
}

\examples{

# load and create example data:
 postNetExample <- postNetParcel::postNetExample
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
                          
# Retrieve the features selected by the model the explain changes in translation efficiency 
 selectedFeatures <- ptn_selectedFeatures(ptn, analysis_type = "lm", comparison = 1) 
 str(selectedFeatures)
 
}
