\docType{methods}
\name{ptn_model}
\alias{ptn_model,postNetData-method}
\title{Accessor for the models feature integration model slots of a \code{postNetData} object
}

\description{
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
  \item{model}{
}
  \item{comparison}{An integer specifying which comparison to select results from. For example, \code{1} would select the first comparison specified in the \code{comparisons} parameter of the \link[postNet]{featureIntegration} function.
}
}

\details{
}

\value{
}

\references{
Sing T, Sander O, Beerenwinkel N, Lengauer T (2005). “ROCR: visualizing classifier performance in R.” Bioinformatics, 21(20), 7881. http://rocr.bioinf.mpi-sb.mpg.de.

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
 omnibusModel <- ptn_model(ptn, analysis_type = "lm", model = "omnibus", comparison = 1) 
 str(omnibusModel)
}
