\docType{methods}
\name{ptn_model}
\alias{ptn_model,postNetData-method}
\title{Accessor for the feature integration model slots of a \code{postNetData} object
}

\description{Retrieve the models resulting from both the stepwise linear regression and random forest implementations of \link[postNet]{featureIntegration} from a \code{postNetData} object.
}

\usage{
\S4method{ptn_model}{postNetData}(ptn, 
                                  analysis_type, 
                                  model, 
                                  comparison)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{analysis_type}{A string specifying the method of analysis used for \link[postNet]{featureIntegration} to retrieve the selected features from. The options are \code{"lm"} to retrieve features selected by stepwise linear regression, or \code{"rf"} for features selected using random forest-based analyses.
}
  \item{model}{A string specifying the model the retrieve. If \code{analysis_type = "lm"}, the options are: \code{"univariateModel"}, \code{"stepwiseModel"}, or \code{"finalModel"}. If \code{analysis_type = "rf"}, the options are: \code{"preModel"}, \code{"borutaModel"}, or \code{"finalModel"}. 
}
  \item{comparison}{An integer specifying which comparison to select results from. For example, \code{1} would select the first comparison specified in the \code{comparisons} parameter of the \link[postNet]{featureIntegration} function.
}
}

\value{If \code{analysis_type = "lm"} and \code{model = "univariateModel"}, and object of class \code{postNetUnivariate} is returned with three slots containing the p-values, FDRs, and proportion of variance explained for each feature passing the significance threshold in univariate models. 

If \code{analysis_type = "lm"} and \code{model = "stepwiseModel"}, and object of class \code{postNetStepWise} is returned with two slots containing the linear models from each step of setpwise regression, and the F-value table storing the F-values for each feature at each step of stepwise regression.

If \code{analysis_type = "lm"} and \code{model = "finalModel"}, and object of class \code{postNetFinalModel} is returned with three slots containing the the total variance explained by the final omnibus model, the final omnibus linear model, and a table summarizing the p-values and variance explained by each feature in the omnibus and adjusted models. 

If \code{analysis_type = "rf"} and \code{model = "preModel"}, and object of class \code{randomForest} is returned for the pre-model using all features. See \link[randomForest]{randomForest} for full description of the output object.

If \code{analysis_type = "rf"} and \code{model = "borutaModel"}, and object of class \code{Boruta} is returned with the results of feature selection. See \link[Boruta]{Boruta} for full description of the output object.

If \code{analysis_type = "rf"} and \code{model = "finalModel"}, and object of class \code{randomForest} is returned for the final model using the selected features. See \link[randomForest]{randomForest} for full description of the output object.
}

\references{
Liaw A, Wiener M (2002). “Classification and Regression by randomForest.” R News, 2(3), 18-22. https://CRAN.R-project.org/doc/Rnews/. 

Witold R. Rudnicki M (2010). “Feature Selection with the Boruta Package.” Journal of Statistical Software, 36(11), 1–13. https://doi.org/10.18637/jss.v036.i11.
}

\seealso{
\link[postNet]{featureIntegration} \cr
\link[Boruta]{Boruta} \cr
\link[randomForest]{randomForest}
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
                          
# Retrieve the omnibus model from stepwise regression analysis
 finalModel <- ptn_model(ptn, analysis_type = "lm", model = "finalModel", comparison = 1) 
 str(finalModel)
}
