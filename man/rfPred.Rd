\name{rfPred}
\alias{rfPred}
\title{Predict post-transcriptional regulation using pre-built random forest models 
}

\description{Aftering identifying features associated with changes in post-transcriptional regulation and building a random forest classification model using the \code{\link[postNet]{featureIntegration}} function, the \code{rfPred} function can then be used to apply this model to predict post-transcriptional regulation in new gene lists or datasets based on the same set of features.  
}

\usage{
rfPred(ptn, 
      comparison, 
      predGeneList, 
      predFeatures, 
      pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{comparison}{An integer specifying which comparison to select results from. For example, \code{1} would select the first comparison specified in the \code{comparisons} parameter of the \link[postNet]{featureIntegration} function.
}
  \item{predGeneList}{A named list of character vectors corresponding to gene sets for which regulation will be predicted. The list must have the same number of gene sets as classes in the predictive model. For example, if the random forest model from \code{\link[postNet]{featureIntegration}} was built to classify genes into activated and suppressed categories, then the input list of genes to make predictions on should correspond to the same categories. Ensure that the gene IDs provided match those in the reference sequence annotations.  
}
  \item{predFeatures}{
}
  \item{pdfName}{Name to be appended to output PDF files.
}
}

\details{
}

\value{
}

\references{
Sing T, Sander O, Beerenwinkel N, Lengauer T (2005). “ROCR: visualizing classifier performance in R.” Bioinformatics, 21(20), 7881. http://rocr.bioinf.mpi-sb.mpg.de.

Liaw A, Wiener M (2002). “Classification and Regression by randomForest.” R News, 2(3), 18-22. https://CRAN.R-project.org/doc/Rnews/. 
}

\seealso{
\code{\link[postNet]{featureIntegraion}} \cr
\code{\link[randomForest]{randomForest}} \cr
\code{\link[ROCR]{ROCR}}
}

\examples{
}
