\name{rfPred}
\alias{rfPred}
\title{Predictions
}

\description{
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
  \item{predGeneList}{
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
