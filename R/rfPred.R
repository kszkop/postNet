rfPred <- function(ptn,
                   comparison,
                   predGeneList,
                   predFeatures,
                   pdfName = NULL) {
  check_ptn(ptn)
  if (is.null(ptn@analysis@featureIntegration$rf)) {
    stop("Please run random forest analysis first")
  }
  if (!check_number(comparison)) {
    stop("Please provide the correct comparison number. You can check them using ptn_check_comparisons(ptn, analysis_type='rf') ")
    if (length(comparison) != 1) {
      stop("There can be only one comparison")
    }
  }
  modelIn <- ptn_model(ptn, analysis_type = "rf", comparison = comparison, model = "finalModel")
  selFeat <- names(ptn_selectedFeatures(ptn, analysis_type = "rf", comparison = comparison))
  #
  check_features(predFeatures)
  if (!all(selFeat %in% names(predFeatures))) {
    missTmp <- setdiff(selFeat, names(predFeatures))
    stop(paste("These features: ", missTmp, "are missing in the predFeature object. Please calculate and add them.", sep = ""))
  }
  if (length(predGeneList) != 2) {
    stop("there can be only 2 entries for predGeneList")
  }

  predFeaturesNames <- names(predFeatures)
  tmpDf <- data.frame(t(plyr::ldply(predFeatures, rbind, .id = NULL)))
  colnames(tmpDf) <- predFeaturesNames
  # Keep only selFeat
  tmpDf <- tmpDf[, colnames(tmpDf) %in% selFeat]

  featIn <- na.omit(tmpDf)
  message(paste(nrow(tmpDf) - nrow(featIn), "genes removed because of NAs", sep = " "))

  listSel <- as.character((unlist(predGeneList)))
  featInSel <- featIn[row.names(featIn) %in% listSel, ]

  featInSel$reg <- NA
  for (i in 1:2) {
    cTmp <- predGeneList[[i]]
    regTmp <- c("A", "B")
    featInSel$reg[row.names(featInSel) %in% cTmp] <- regTmp[i]
  }

  featInSel <- featInSel[!is.na(featInSel$reg), ]
  featInSel$reg <- as.factor(featInSel$reg)
  #
  predValidc <- stats::predict(modelIn, featInSel, type = "class")
  predValid <- stats::predict(modelIn, featInSel, type = "prob")
  #
  perf <- ROCR::prediction(predValid[, 2], as.numeric(featInSel$reg))
  # 1. Area under curve
  # auc = ROCR::performance(perf, "auc")
  # 2. True Positive and Negative Rate
  predOut <- ROCR::performance(perf, "tpr", "fpr")
  # 3. Plot the ROC curve
  pdf(ifelse(is.null(pdfName), "rocr.pdf", paste(pdfName, "_rocr.pdf", sep = "")), width = 8, height = 8, useDingbats = F)
  plot(predOut, main = paste("ROC Curve for Random Forest \n Accuracy: ", round(as.numeric(caret::confusionMatrix(predValidc, featInSel$reg)[[3]][1]), 3), sep = ""), col = "firebrick1", lwd = 3, xlab = "", ylab = "", )
  abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")

  mtext(side = 1, line = 4, "False positive rate", col = "black", font = 2, cex = 1.2)
  mtext(side = 2, line = 3, "True positive rate", col = "black", font = 2, cex = 1.2)
  text(0.8, 0.2, font = 2, cex = 1.7, paste("Sensitivity: ", round(caret::confusionMatrix(predValidc, featInSel$reg)[[4]][1], 2), sep = ""))
  text(0.8, 0.1, font = 2, cex = 1.7, paste("Specificity: ", round(caret::confusionMatrix(predValidc, featInSel$reg)[[4]][2], 2), sep = ""))
  dev.off()

  # ptn@analysis@featureIntegration@rf@prediction <- perf
  return(ptn)
}
