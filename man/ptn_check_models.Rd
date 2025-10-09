\name{ptn_check_models}
\alias{ptn_check_models}
\title{Accessor for the feature integration model comparisons in a \code{postNetData} object
}

\description{Retrieve the comparisons used for the models resulting from both the stepwise linear regression and random forest implementations of \link[postNet]{featureIntegration} from a \code{postNetData} object.
}

\usage{
ptn_check_models(ptn, analysis_type)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{analysis_type}{A string specifying the method of analysis used for \link[postNet]{featureIntegration} to retrieve the selected features from. The options are \code{"lm"} to retrieve features selected by stepwise linear regression, or \code{"rf"} for features selected using random forest-based analyses.
}
}

\value{Returns a character string specifying the comparisons that were used in modeling. For example "translationUp_c1_translationDown_c1" indicates the gene sets that were compared, and from which anota2seq contrast they were selected from (e.g. c1 = contrast 1).
}

\seealso{
\link[postNet]{featureIntegration} \cr
\link[postNet]{ptn_model}
}

\examples{# load and create example data:
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
                          
# Retrieve the contrasts that were used in modelling
  ptn_check_models(ptn, analysis_type = "lm") 
}
