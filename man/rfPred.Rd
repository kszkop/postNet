\name{rfPred}
\alias{rfPred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rfPred(ptn, comparison, predGeneList, predFeatures, pdfName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{comparison}{
%%     ~~Describe \code{comparison} here~~
}
  \item{predGeneList}{
%%     ~~Describe \code{predGeneList} here~~
}
  \item{predFeatures}{
%%     ~~Describe \code{predFeatures} here~~
}
  \item{pdfName}{
%%     ~~Describe \code{pdfName} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (ptn, comparison, predGeneList, predFeatures, pdfName = NULL) 
{
    check_ptn(ptn)
    if (is.null(slot(ptn@analysis@featureIntegration, "rf"))) {
        stop(paste("Please run ", analysis_type, "analysis first", 
            sep = ""))
    }
    if (!check_number(comparison)) {
        stop("please provide correct comparison number. You can check them using ptn_check_comparisons(ptn, analysis_type='rf') ")
        if (length(comparison) != 1) {
            stop("there can be only one comparison")
        }
    }
    modelIn <- ptn_model(ptn, analysis_type = "rf", comparison = comparison, 
        model = "finalModel")
    selFeat <- names(ptn_selectedFeatures(ptn, analysis_type = "rf", 
        comparison = comparison))
    check_features(predFeatures)
    if (!all(selFeat \%in\% names(predFeatures))) {
        missTmp <- setdiff(selFeat, names(predFeatures))
        stop(paste("These features: ", missTmp, "are missing in the predFeature object. Please calculate and add them.", 
            sep = ""))
    }
    if (length(predGeneList) != 2) {
        stop("there can be only 2 entries for predGeneList")
    }
    predFeaturesNames <- names(predFeatures)
    tmpDf <- data.frame(t(plyr::ldply(predFeatures, rbind, .id = NULL)))
    colnames(tmpDf) <- predFeaturesNames
    tmpDf <- tmpDf[, colnames(tmpDf) \%in\% selFeat]
    featIn <- na.omit(tmpDf)
    message(paste(nrow(tmpDf) - nrow(featIn), "genes removed because of NAs", 
        sep = " "))
    listSel <- as.character((unlist(predGeneList)))
    featInSel <- featIn[row.names(featIn) \%in\% listSel, ]
    featInSel$reg <- NA
    for (i in 1:2) {
        cTmp <- predGeneList[[i]]
        regTmp <- c("A", "B")
        featInSel$reg[row.names(featInSel) \%in\% cTmp] <- regTmp[i]
    }
    featInSel <- featInSel[!is.na(featInSel$reg), ]
    featInSel$reg <- as.factor(featInSel$reg)
    predValidc <- stats::predict(modelIn, featInSel, type = "class")
    predValid <- stats::predict(modelIn, featInSel, type = "prob")
    perf = ROCR::prediction(predValid[, 2], as.numeric(featInSel$reg))
    predOut = ROCR::performance(perf, "tpr", "fpr")
    pdf("eif4e_rocr.pdf", width = 8, height = 8, useDingbats = F)
    plot(predOut, main = paste("ROC Curve for Random Forest \n Accuracy: ", 
        round(as.numeric(caret::confusionMatrix(predValidc, featInSel$reg)[[3]][1]), 
            3), sep = ""), col = "firebrick1", lwd = 3, xlab = "", 
        ylab = "", )
    abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")
    mtext(side = 1, line = 4, "False positive rate", col = "black", 
        font = 2, cex = 1.2)
    mtext(side = 2, line = 3, "True positive rate", col = "black", 
        font = 2, cex = 1.2)
    text(0.8, 0.2, font = 2, cex = 1.7, paste("Sensitivity: ", 
        round(caret::confusionMatrix(predValidc, featInSel$reg)[[4]][1], 
            2), sep = ""))
    text(0.8, 0.1, font = 2, cex = 1.7, paste("Specificity: ", 
        round(caret::confusionMatrix(predValidc, featInSel$reg)[[4]][2], 
            2), sep = ""))
    dev.off()
    return(ptn)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
