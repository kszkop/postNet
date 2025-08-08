\name{rfPred}
\alias{rfPred}
\title{Predict post-transcriptional regulation using pre-built random forest models 
}

\description{Aftering identifying features associated with changes in post-transcriptional regulation and training a random forest classification model using the \code{\link[postNet]{featureIntegration}} function, the \code{rfPred} function can then be used to apply this model to predict post-transcriptional regulation in new gene lists or datasets based on the same set of features. In this way, it is possible to evaluate how well a set of features identified as being associated with post-transcriptional regulation in one context can explain regulation observed in another context.   
}

\usage{
rfPred(ptn, 
      comparison, 
      predGeneList, 
      predFeatures, 
      pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object after running the random forest implementation of \code{\link[postNet]{featureIntegration}},  containing the trained final random forest model.
}
  \item{comparison}{An integer specifying which comparison to select the model from. For example, \code{1} would select the first comparison specified in the \code{comparisons} parameter of the \link[postNet]{featureIntegration} function.
}
  \item{predGeneList}{A named list of character vectors corresponding to the new gene lists for which regulation will be predicted. The list must have the same number of gene lists as classes in the predictive model trained on the original data. For example, if the random forest model from \code{\link[postNet]{featureIntegration}} was trained to classify genes into activated and suppressed categories, then the input list of genes to make predictions on should correspond to the same categories in the new dataset. Ensure that the gene IDs provided match those in the reference sequence annotations.  
}
  \item{predFeatures}{A named list of numeric vectors corresponding to features that will be used with the model trained on the original data to predict regulatory effects in the new data. Each vector element in the list must have names corresponding to gene IDs. These features must be enumerated in the new gene lists or dataset that is being predicted on, and must include the same features used in the final random forest model trained on the original dataset.   
}
  \item{pdfName}{Name to be appended to output PDF files.
}
}

\details{
}

\value{The output of \code{rfPred} is a PDF file with the Receiver Operating Characteristic (ROC) curve illustrating the performance of the trained model in predicting in the new dataset.
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
