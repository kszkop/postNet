\name{featureIntegration}
\alias{featureIntegration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
featureIntegration(ptn, features, lmfeatGroup = NULL, lmfeatGroupColour = NULL, analysis_type, regOnly = TRUE, allFeat = FALSE, useCorel = TRUE, covarFilt = 20, NetModelSel = "omnibus", comparisons = NULL, pdfName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{features}{
%%     ~~Describe \code{features} here~~
}
  \item{lmfeatGroup}{
%%     ~~Describe \code{lmfeatGroup} here~~
}
  \item{lmfeatGroupColour}{
%%     ~~Describe \code{lmfeatGroupColour} here~~
}
  \item{analysis_type}{
%%     ~~Describe \code{analysis_type} here~~
}
  \item{regOnly}{
%%     ~~Describe \code{regOnly} here~~
}
  \item{allFeat}{
%%     ~~Describe \code{allFeat} here~~
}
  \item{useCorel}{
%%     ~~Describe \code{useCorel} here~~
}
  \item{covarFilt}{
%%     ~~Describe \code{covarFilt} here~~
}
  \item{NetModelSel}{
%%     ~~Describe \code{NetModelSel} here~~
}
  \item{comparisons}{
%%     ~~Describe \code{comparisons} here~~
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
function (ptn, features, lmfeatGroup = NULL, lmfeatGroupColour = NULL, 
    analysis_type, regOnly = TRUE, allFeat = FALSE, useCorel = TRUE, 
    covarFilt = 20, NetModelSel = "omnibus", comparisons = NULL, 
    pdfName = NULL) 
{
    check_ptn(ptn)
    check_features(features)
    if (!is.null(comparisons)) {
        if (!check_comparisons(comparisons)) {
            stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
        }
        if (length(which(unique(unlist(comparisons)) == 0)) > 
            0 && is.null(ptn_background(ptn))) {
            stop(" 0 is always a background, but no background provided")
        }
    }
    check_analysis_type(analysis_type)
    if (!check_logical(regOnly)) {
        stop("'regOnly' can only be TRUE or FALSE")
    }
    if (analysis_type == "lm") {
        if (!check_logical(allFeat)) {
            stop("'allFeat' can only be TRUE or FALSE")
        }
        if (!check_logical(useCorel)) {
            stop("'useCorel' can only be TRUE or FALSE")
        }
        if (!check_number(covarFilt)) {
            stop("'covarFilt' can only be a numerical value")
        }
        if (!is_valid_NetModelSel(NetModelSel)) {
            stop("'NetModelSel' has to be not null and one one of the: omnibus, adjusted or univariate")
        }
    }
    if (analysis_type == "rf") {
        if (is.null(comparisons)) {
            stop("Please provide desired comparisons")
        }
    }
    ptn <- prepFeatures(ptn, features)
    dataTmp <- ptn_features(ptn)
    effTmp <- ptn_effect(ptn)
    colnames(dataTmp) <- c(paste("a", seq(1, ncol(dataTmp), 1), 
        sep = ""))
    dataTmp$effM <- effTmp[match(row.names(dataTmp), names(effTmp))]
    namesDf <- data.frame(originalNames = colnames(ptn_features(ptn))[1:ncol(dataTmp) - 
        1], newNames = colnames(dataTmp)[1:ncol(dataTmp) - 1], 
        stringsAsFactors = F)
    resOut <- resQuant(qvec = ptn_effect(ptn), ptn = ptn)
    fiOut <- new("postNetFeatureIntegration", lm = NULL, rf = NULL, 
        featureMap = NULL)
    if (analysis_type == "lm") {
        if (!is.null(lmfeatGroup)) {
            checklmfeatGroup(lmfeatGroup, ncol(dataTmp) - 1)
            names(lmfeatGroup) <- colnames(dataTmp)[1:ncol(dataTmp) - 
                1]
            checklmfeatGroupColour(lmfeatGroupColour, lmfeatGroup)
            names(lmfeatGroupColour) <- lmfeatGroup
        }
        if (isTRUE(regOnly)) {
            for (i in 1:length(comparisons)) {
                coloursTmp <- ptn_colours(ptn)
                if (names(resOut)[1] == "background") {
                  compTmp <- comparisons[[i]] + 1
                  coloursTmp <- c("grey75", coloursTmp)[compTmp]
                }
                else {
                  compTmp <- comparisons[[i]]
                  coloursTmp <- coloursTmp[compTmp]
                }
                listSel <- c(names(resOut[[compTmp[1]]]), names(resOut[[compTmp[2]]]))
                dataTmpSel <- dataTmp[row.names(dataTmp) \%in\% 
                  listSel, ]
                lmOut <- runLM(dataIn = dataTmpSel, namesDf = namesDf, 
                  allFeat = allFeat, useCorel = useCorel, covarFilt = covarFilt, 
                  nameOut = pdfName, NetModelSel = NetModelSel, 
                  coloursIn = coloursTmp, lmfeatGroup = lmfeatGroup, 
                  lmfeatGroupColour = lmfeatGroupColour)
                fiOut@lm[[paste(names(resOut)[compTmp], collapse = "_")]] <- lmOut
                bestSel <- names(lmOut@selectedFeatures)
                for (feat in bestSel) {
                  featTmp <- namesDf[namesDf$originalNames == 
                    feat, ]$newNames
                  set <- dataTmpSel[, colnames(dataTmpSel) \%in\% 
                    c(featTmp, "effM")]
                  set1 <- names(resOut[[compTmp[1]]])
                  setSel1 <- set[row.names(set) \%in\% set1, ]
                  set2 <- names(resOut[[compTmp[2]]])
                  setSel2 <- set[row.names(set) \%in\% set2, ]
                  plotScatterInd(set1 = setSel1, set2 = setSel2, 
                    orgName = feat, coloursIn = coloursTmp, nameOut = pdfName)
                }
            }
        }
        else {
            coloursTmp <- c("salmon", "skyblue")
            lmOut <- runLM(dataIn = dataTmp, namesDf = namesDf, 
                allFeat = allFeat, useCorel = useCorel, covarFilt = covarFilt, 
                nameOut = pdfName, NetModelSel = NetModelSel, 
                coloursIn = coloursTmp, lmfeatGroup = lmfeatGroup, 
                lmfeatGroupColour = lmfeatGroupColour)
            fiOut@lm[["allData"]] <- lmOut
            bestSel <- names(lmOut@selectedFeatures)
            for (feat in bestSel) {
                featTmp <- namesDf[namesDf$originalNames == feat, 
                  ]$newNames
                set <- dataTmp[, colnames(dataTmp) \%in\% c(featTmp, 
                  "effM")]
                plotScatterInd(set1 = set, set2 = NULL, orgName = feat, 
                  coloursIn = "grey75", nameOut = pdfName)
            }
        }
    }
    else if (analysis_type == "rf") {
        dataTmpReg <- dataTmp[, colnames(dataTmp) != "effM"]
        colnames(dataTmpReg) <- namesDf$originalNames[match(colnames(dataTmpReg), 
            namesDf$newNames)]
        for (i in 1:length(comparisons)) {
            coloursTmp <- ptn_colours(ptn)
            if (names(resOut)[1] == "background") {
                compTmp <- comparisons[[i]] + 1
                coloursTmp <- c("grey75", coloursTmp)[compTmp]
            }
            else {
                compTmp <- comparisons[[i]]
                coloursTmp <- coloursTmp[compTmp]
            }
            dataTmpSel <- dataTmpReg
            dataTmpSel$reg <- NA
            for (j in 1:2) {
                if (length(compTmp) != 2) {
                  stop("There is something wrong with comparisons")
                }
                cTmp <- names(resOut[[compTmp[j]]])
                regTmp <- c("A", "B")
                dataTmpSel$reg[row.names(dataTmpSel) \%in\% cTmp] <- regTmp[j]
            }
            dataTmpSel <- dataTmpSel[!is.na(dataTmpSel$reg), 
                ]
            dataTmpSel$reg <- as.factor(dataTmpSel$reg)
            train <- sample(nrow(dataTmpSel), 0.7 * nrow(dataTmpSel), 
                replace = FALSE)
            TrainSet <- dataTmpSel[train, ]
            TrainSet$reg <- as.factor(TrainSet$reg)
            ValidSet <- dataTmpSel[-train, ]
            ValidSet$reg <- as.factor(ValidSet$reg)
            model1 <- randomForest::randomForest(reg ~ ., data = TrainSet, 
                importance = TRUE, ntree = 500)
            model1Imp <- Boruta::Boruta(reg ~ ., data = TrainSet, 
                doTrace = 0, maxRuns = 500, pValue = 0.001)
            featComf <- row.names(Boruta::attStats(model1Imp))[which(as.character(Boruta::attStats(model1Imp)[, 
                6]) == "Confirmed")]
            pdf(paste("featureImportance.pdf", sep = "_"), width = 8, 
                height = 8, useDingbats = F)
            par(mar = c(10, 5, 3, 3), bty = "l", font = 2, font.axis = 2, 
                font.lab = 2, cex.axis = 1.3, cex.main = 1.7, 
                cex.lab = 1)
            plot(model1Imp, las = 2, xlab = "", ylab = "", yaxt = "n", 
                xaxt = "n", pch = 20)
            mtext(side = 1, line = 9, "Features", col = "black", 
                font = 2, cex = 1.2)
            mtext(side = 2, line = 3, "Importance (Z-score)", 
                col = "black", font = 2, cex = 1.2)
            axis(side = 2, seq(0, roundNice(max(Boruta::attStats(model1Imp)[, 
                4]), direction = "up"), 10), font = 2, lwd = 2, 
                las = 2, cex = 0.75)
            tmp <- lapply(1:ncol(model1Imp$ImpHistory), function(i) model1Imp$ImpHistory[is.finite(model1Imp$ImpHistory[, 
                i]), i])
            names(tmp) <- colnames(model1Imp$ImpHistory)
            tmpNames <- names(sort(sapply(tmp, median)))
            addNames <- c("shadowMin", "shadowMax", "shadowMean")
            coloursN <- rep("black", length(tmpNames))
            coloursN[tmpNames \%in\% addNames] <- "firebrick1"
            axis(side = 1, at = 1:length(tmpNames), labels = F, 
                font = 2, lwd = 2, las = 2, cex.axis = 0.5, tck = -0.005)
            text(1:length(tmpNames), par("usr")[3] - 1.05, labels = tmpNames, 
                col = coloursN, srt = 45, adj = 1, cex = 0.55, 
                xpd = NA)
            dev.off()
            TrainSet <- TrainSet[, colnames(TrainSet) \%in\% c(featComf, 
                "reg")]
            ValidSet <- ValidSet[, colnames(ValidSet) \%in\% c(featComf, 
                "reg")]
            model2 <- randomForest::randomForest(reg ~ ., data = TrainSet, 
                importance = TRUE, ntree = 500)
            varImpIn <- sort(randomForest::importance(model2)[, 
                3], decreasing = T)
            pdf(paste("FinalModel.pdf", sep = "_"), width = 16, 
                height = 8, useDingbats = F)
            par(mfrow = c(1, 2), mar = c(9, 5, 10, 4), bty = "l", 
                font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, 
                cex.main = 1.7, cex.lab = 1)
            colDot <- rep("black", length(randomForest::importance(model2)[, 
                3]))
            dotchart(sort(randomForest::importance(model2)[, 
                3], decreasing = F), cex = 0.75, col = colDot, 
                labels = names(sort(varImpIn, decreasing = F)), 
                xlab = "", xaxt = "n", frame.plot = FALSE, pch = 16)
            axis(side = 1, seq(0, roundNice(max(varImpIn), direction = "up"), 
                5), font = 2, lwd = 2)
            mtext(side = 1, line = 4, "Feature Importance \n (Mean Decrease Accuracy)", 
                col = "black", font = 2, cex = 1.2)
            predValidc <- stats::predict(model2, ValidSet, type = "class")
            predValid <- stats::predict(model2, ValidSet, type = "prob")
            perf <- ROCR::prediction(predValid[, 2], as.numeric(ValidSet$reg))
            auc <- ROCR::performance(perf, "auc")
            predOut <- ROCR::performance(perf, "tpr", "fpr")
            plot(predOut, main = paste("ROC Curve for Random Forest \n AUC: ", 
                round(auc@y.values[[1]], 3), sep = ""), col = "firebrick1", 
                lwd = 3, xlab = "", ylab = "", )
            abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")
            mtext(side = 1, line = 4, "False positive rate", 
                col = "black", font = 2, cex = 1.2)
            mtext(side = 2, line = 3, "True positive rate", col = "black", 
                font = 2, cex = 1.2)
            text(0.8, 0.2, font = 2, cex = 1.7, paste("Sensitivity: ", 
                round(caret::confusionMatrix(predValidc, ValidSet$reg)[[4]][1], 
                  2), sep = ""))
            text(0.8, 0.1, font = 2, cex = 1.7, paste("Specificity: ", 
                round(caret::confusionMatrix(predValidc, ValidSet$reg)[[4]][2], 
                  2), sep = ""))
            dev.off()
            rfOut <- new("postNetFeatureIntegration_rf", preModel = model1, 
                borutaModel = model1Imp, finalModel = model2, 
                selectedFeatures = varImpIn)
            fiOut@rf[[paste(names(resOut)[compTmp], collapse = "_")]] <- rfOut
            bestSel <- names(rfOut@selectedFeatures)
            for (feat in bestSel) {
                featTmp <- namesDf[namesDf$originalNames == feat, 
                  ]$newNames
                set <- ptn_effect(ptn)
                set <- dataTmp[, colnames(dataTmp) \%in\% c(featTmp, 
                  "effM")]
                set1 <- names(resOut[[compTmp[1]]])
                setSel1 <- set[row.names(set) \%in\% set1, ]
                set2 <- names(resOut[[compTmp[2]]])
                setSel2 <- set[row.names(set) \%in\% set2, ]
                plotScatterInd(set1 = setSel1, set2 = setSel2, 
                  orgName = feat, coloursIn = coloursTmp, nameOut = pdfName)
            }
        }
    }
    else {
        stop("Please provide correct type: lm for linear regression or rf for random forest")
    }
    ptn@analysis@featureIntegration <- fiOut
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
