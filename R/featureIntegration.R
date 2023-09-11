featureIntegration <- function(features,
                               analysis_type,
                               rmCat = NULL,
                               regOnly = TRUE,
                               allFeat = TRUE,
                               useCorel=TRUE,
                               covarFilt = 20,
                               NetModelSel = "Omnibus",
                               ads = NULL,
                               regulation = NULL,
                               contrast = NULL,
                               geneList = NULL,
                               geneListcolours = NULL,
                               customBg = NULL,
                               comparisons = NULL,
                               regulationGen = NULL,
                               contrastSel = NULL,
                               effectMeasure = NULL,
                               outDir = NULL,
                               pdfName = NULL) {
  #
  if (!is.null(outDir)){
    dirTmp <- outDir
  } else {
    dirTmp <- paste('featureIntegration', format(Sys.time(), "%Y%m%e_%X"),sep='_')
  }
  nameTmp <- ifelse(is.null(pdfName), paste(ifelse(is.null(regulationGen), "custom", regulationGen), sep = "_"), paste(pdfName, ifelse(is.null(regulationGen), "custom", regulationGen), sep = "_"))
  #
  if (!is.null(effectMeasure)) {
    effM <- effectMeasure
  } else {
    if(regulationGen=='mRNAAbundance'){
      regTmp <- 'totalmRNA'
    } else {
      regTmp <- regulationGen
    }
    scOut <- anota2seq::anota2seqGetOutput(ads, output = "singleDf", selContrast = contrastSel, getRVM = TRUE)
    effM <- scOut[, grepl(paste(regTmp, "apvEff", sep = "."), colnames(scOut))]
    names(effM) <- scOut$identifier
  }
  #
  featureNames <- names(features)
  featuresTmp <- append(features, list(effM))
  featuresTmp <- unname(featuresTmp)
  
  tmpDf <- data.frame(t(plyr::ldply(featuresTmp, rbind)))
  tmpDf <- na.omit(tmpDf)
  
  dat <- tmpDf
  colnames(dat) <- c(paste("a", seq(1, length(featureNames), 1), sep = ""), "effM")
  
  namesDf <- data.frame(originalNames = featureNames, newNames = colnames(dat)[-length(colnames(dat))], stringsAsFactors = F)
  dataOrg <- dat

  resOut <- resSel(vIn = effM, ads = ads, regulation = regulation, contrast = contrast, customBg = customBg, geneList = geneList)

  ######
  if (analysis_type == 'lm'){
    if (isTRUE(regOnly)){
      for (i in 1:length(comparisons)) {
        if (!is.null(ads) | !is.null(customBg)) {
          compTmp <- comparisons[[i]] + 1
        } else {
          compTmp <- comparisons[[i]]
        }
        #
        listSel <- c(names(resOut[[compTmp[1]]]), names(resOut[[compTmp[2]]]))
        dat <- dat[row.names(dat) %in% listSel, ]
        #
        dirTmp <- paste(dirTmp, "lm", paste(compTmp[1], "vs", compTmp[2], sep = "_"),sep='_')
        dir.create(dirTmp)
        nameOut <- paste(dirTmp,nameTmp, sep='/')
        #
        runLM(dataIn = dat, namesDf = namesDf, allFeat = allFeat, useCorel = useCorel, nameOut = nameOut, NetModelSel = NetModelSel)
      }
    } else {
      dirTmp <- paste(dirTmp, 'lm',sep='_')
      dir.create(dirTmp)
      nameOut <- paste(dirTmp,nameTmp, sep='/')
      #
      runLM(dataIn = dat, namesDf = namesDf, allFeat = allFeat, useCorel = useCorel, nameOut = nameOut, NetModelSel = NetModelSel)
    }
  } else if (analysis_type == "rf") {
    #
    dirTmp <- paste(dirTmp, "randomforest",sep='_')
    dir.create(dirTmp)
    nameOut <- paste(dirTmp,nameTmp, sep='/')
    #
    if (!is.null(ads) | !is.null(customBg)) {
      catTmp <- rmCat + 1
    } else {
      catTmp <- rmCat
    }
    dat$reg <- NA
    for(i in catTmp) {
      regTmp <- names(resOut[[i]])
      dat$reg[row.names(dat) %in% regTmp] <- i
    }
    dat <- dat[!is.na(dat$reg),]
    #
    dat <- dat[, colnames(dat) != "effM"]
    # Split to Train and Valid ( 70:30)
    train <- sample(nrow(dat), 0.7 * nrow(dat), replace = FALSE)
    TrainSet <- dat[train, ]
    TrainSet$reg <- as.factor(TrainSet$reg)
    ValidSet <- dat[-train, ]
    ValidSet$reg <- as.factor(ValidSet$reg)
    # run model on training set
    model1 <- randomForest::randomForest(reg ~ ., data = TrainSet, importance = TRUE, ntree = 500)
    
    # Plot importance
    # also apply to find relevant features
    model1Imp <- Boruta::Boruta(reg ~ ., data = TrainSet, doTrace = 0, maxRuns = 500, pValue = 0.001)
    # selecct important once
    featComf <- row.names(Boruta::attStats(model1Imp))[which(as.character(Boruta::attStats(model1Imp)[, 6]) == "Confirmed")]
    #
    pdf(paste(nameOut, "featureImportance.pdf", sep = "_"), width = 8, height = 8, useDingbats = F)
    par(mar = c(10, 5, 3, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
    plot(model1Imp, las = 2, xlab = "", ylab = "", yaxt = "n", xaxt = "n", pch = 20)
    mtext(side = 1, line = 9, "Features", col = "black", font = 2, cex = 1.2)
    mtext(side = 2, line = 3, "Importance (Z-score)", col = "black", font = 2, cex = 1.2)
    axis(side = 2, seq(0, roundUpNice(max(Boruta::attStats(model1Imp)[, 4])), 10), font = 2, lwd = 2, las = 2, cex = 0.75)
    #
    tmp <- lapply(1:ncol(model1Imp$ImpHistory), function(i) model1Imp$ImpHistory[is.finite(model1Imp$ImpHistory[, i]), i])
    names(tmp) <- colnames(model1Imp$ImpHistory)
    tmpNames <- names(sort(sapply(tmp, median)))
    addNames <- c("shadowMin", "shadowMax", "shadowMean")
    
    tmpNames <- c(namesDf$originalNames, addNames)[match(tmpNames, c(namesDf$newNames, addNames))]
    coloursN <- rep("black", length(tmpNames))
    coloursN[tmpNames %in% addNames] <- "firebrick1"
    axis(side = 1, at = 1:length(tmpNames), labels = F, font = 2, lwd = 2, las = 2, cex.axis = 0.5)
    text(1:length(tmpNames), par("usr")[3] - 1.05, labels = tmpNames, col = coloursN, srt = 45, adj = 1, cex = 0.55, xpd = NA)
    dev.off()
      
    # rerun model with only confirmed ones
    TrainSet <- TrainSet[, colnames(TrainSet) %in% c(featComf, "reg")]
    ValidSet <- ValidSet[, colnames(ValidSet) %in% c(featComf, "reg")]
      
    # run model on training set
    model2 <- randomForest::randomForest(reg ~ ., data = TrainSet, importance = TRUE, ntree = 500)
    #
    varImpIn <- sort(randomForest::importance(model1)[, 3], decreasing = T)
    #
    pdf(paste(nameOut,"FinalModel.pdf", sep = "_"), width = 16, height = 8, useDingbats = F)
    par(mfrow = c(1, 2), mar = c(9, 5, 10, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
    colDot <- rep("black", length(randomForest::importance(model1)[, 3]))
    colDot[which(names(sort(randomForest::importance(model1)[, 3], decreasing = F)) %in% featComf)] <- "#B0F2BC"
    dotchart(sort(randomForest::importance(model1)[, 3], decreasing = F), cex = 0.75, col = colDot, labels = namesDf$originalNames[match(names(sort(randomForest::importance(model1)[, 3], decreasing = F)), namesDf$newNames)], xlab = "", xaxt = "n", frame.plot = FALSE, pch = 16)
        
    axis(side = 1, seq(0, roundUpNice(max(varImpIn)), 5), font = 2, lwd = 2)
    mtext(side = 1, line = 4, "Feature Importance \n (Mean Decrease Accuracy)", col = "black", font = 2, cex = 1.2)
        
    # conf <- model1$confusion[,-ncol(model1$confusion)]
    # oob <- (1 - (sum(diag(conf))/sum(conf)))*100
    # mtext(side=3, line=3, paste('OOB estimate of  error rate: ',round(oob,2),sep=''), col="black", font=2, cex=1.7)
        
    # Mean Decrease Accuracy - How much the model accuracy decreases if we drop that variable.
    # Mean Decrease Gini - Measure of variable importance based on the Gini impurity index used for the calculation of splits in trees.
    predValidc <- stats::predict(model2, ValidSet, type = "class")
    #
    predValid <- stats::predict(model2, ValidSet, type = "prob")
      #
    perf <- ROCR::prediction(predValid[, 2], as.numeric(ValidSet$reg))
    # 1. Area under curve
    auc <- ROCR::performance(perf, "auc")
    # 2. True Positive and Negative Rate
    predOut <- ROCR::performance(perf, "tpr", "fpr")
    # 3. Plot the ROC curve
    plot(predOut, main = paste("ROC Curve for Random Forest \n AUC: ", round(auc@y.values[[1]], 3), sep = ""), col = "firebrick1", lwd = 3, xlab = "", ylab = "", )
    abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")
        
    mtext(side = 1, line = 4, "False positive rate", col = "black", font = 2, cex = 1.2)
    mtext(side = 2, line = 3, "True positive rate", col = "black", font = 2, cex = 1.2)
    text(0.8, 0.2, font = 2, cex = 1.7, paste("Sensitivity: ", round(caret::confusionMatrix(predValidc, ValidSet$reg)[[4]][1], 2), sep = ""))
    text(0.8, 0.1, font = 2, cex = 1.7, paste("Specificity: ", round(caret::confusionMatrix(predValidc, ValidSet$reg)[[4]][2], 2), sep = ""))
    dev.off()
  } else {
    stop("Please provide correct type: lm for linear regression or rf for random forest")
  }
  #
  # Prepare plotting
  coloursOut <- coloursSel(ads = ads, regulationGen=regulationGen, geneList = geneList, geneListcolours = geneListcolours, customBg = customBg)
  #
  if (analysis_type == "rf") {
    bestSel <- featComf
  }
  for (feat in bestSel) {
    #
    featTmp <- namesDf[namesDf$newNames == feat, ]$originalNames
    featTmp <- gsub(" ", "_", featTmp)
    #
    if (isTRUE(regOnly) | !is.null(geneList)) {
      #
      set <- dataOrg[row.names(dataOrg) %in% listSel, colnames(dataOrg) == feat]
      set_effm <- dataOrg$effM[row.names(dataOrg) %in% listSel]
      #
      if (is.null(geneList)) {
        set1 <- dataOrg[row.names(dataOrg) %in% as.character(resTmp[grepl(regulationGen, names(resTmp))][[1]]), colnames(dataOrg) == feat]
        set2 <- dataOrg[row.names(dataOrg) %in% as.character(resTmp[grepl(regulationGen, names(resTmp))][[2]]), colnames(dataOrg) == feat]
        set_effm1 <- dataOrg$effM[row.names(dataOrg) %in% as.character(resTmp[grepl(regulationGen, names(resTmp))][[1]])]
        set_effm2 <- dataOrg$effM[row.names(dataOrg) %in% as.character(resTmp[grepl(regulationGen, names(resTmp))][[2]])]
        col1 <- coloursOut[2]
        col2 <- coloursOut[3]
      } else {
        set1 <- dataOrg[row.names(dataOrg) %in% as.character(unlist(resTmp[c(1)])), colnames(dataOrg) == feat]
        set2 <- dataOrg[row.names(dataOrg) %in% as.character(unlist(resTmp[c(2)])), colnames(dataOrg) == feat]
        set_effm1 <- dataOrg$effM[row.names(dataOrg) %in% as.character(unlist(resTmp[c(1)]))]
        set_effm2 <- dataOrg$effM[row.names(dataOrg) %in% as.character(unlist(resTmp[c(2)]))]
        col1 <- geneListcolours[1]
        col2 <- geneListcolours[2]
      }
      
      #
    } else {
      #
      set <- dataOrg[, colnames(dataOrg) == feat]
      set_effm <- dataOrg$effM
      #
    }
    pdf(paste(nameOut, featTmp, "individually.pdf", sep = "_"), width = 8, height = 8, useDingbats = F)
    par(mar = c(9, 5, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.3)
    
    #
    xlim_max <- ifelse(max(set) >= 0, roundUpNice(abs(max(set))), -roundUpNice(abs(max(set))))
    if (xlim_max < max(set)) {
      xlim_max <- ceiling(max(set))
    }
    xlim_min <- ifelse(min(set) >= 0, roundUpNice(abs(min(set))), -roundUpNice(abs(min(set))))
    if (xlim_min > min(set)) {
      xlim_min <- floor(min(set))
    }
    ylim_max <- roundUpNice(abs(max(set_effm)))
    #
    plot(set, set_effm, col = "#8A8683", pch = 16, cex = 1, ylim = c(-ylim_max, ylim_max), xlab = "", ylab = "", lwd = 1, bty = "n", xaxt = "n", yaxt = "n", font = 2, xlim = c(xlim_min, xlim_max))
    
    if (isTRUE(regOnly) | !is.null(geneList)) {
      points(set1, set_effm1, pch = 16, col = col1)
      points(set2, set_effm2, pch = 16, col = col2)
    }
    
    #
    mtext(side = 2, line = 3, ifelse(is.null(regulationGen), "custom", regulationGen), col = "black", font = 2, cex = 1.7)
    axis(side = 2, seq(-ylim_max, ylim_max, 2), font = 2, las = 2, lwd = 2)
    
    mtext(side = 1, line = 4, featTmp, col = "black", font = 2, cex = 1.7, at = (xlim_min + xlim_max) / 2)
    axis(side = 1, seq(xlim_min, xlim_max, ifelse((xlim_max - xlim_min) / 5 >= 0, roundUpNice((xlim_max - xlim_min) / 5), -roundUpNice(abs((xlim_max - xlim_min) / 5)))), font = 2, lwd = 2)
    #
    if (length(unique(set)) > 2 & IQR(set) > 0 & length(unique(set)) > 3) {
      f1 <- predict(smooth.spline(set_effm ~ set))
      lines(f1$x[which(f1$x > xlim_min & f1$x < xlim_max)], f1$y[which(f1$x > xlim_min & f1$x < xlim_max)], col = "#AFBADC", lwd = 4, lend = 2)
      lines(f1$x[which(f1$x > xlim_min & f1$x < xlim_max)], f1$y[which(f1$x > xlim_min & f1$x < xlim_max)], col = "black", lwd = 1, lend = 2, lty = 3)
    }
    #
    if (!is.na(as.numeric(coefficients(lm(set_effm ~ set))[2]))) {
      plotrix::ablineclip(lm(set_effm ~ set), col = "#AFBADC", lwd = 4, x1 = xlim_min, x2 = xlim_max)
      plotrix::ablineclip(lm(set_effm ~ set), col = "black", lwd = 1, x1 = xlim_min, x2 = xlim_max)
      
      text((xlim_min + xlim_max) / 2, ylim_max, paste("pvalue ", format(as.numeric(cor.test(set, set_effm)[3]), scientific = T, digits = 3), ", r=", round(as.numeric(cor.test(set, set_effm)[4]), 3), sep = ""), bty = "n", col = "black", cex = 1.25, font = 2)
    }
    dev.off()
  }
  
  # heatmap select only important/omnibus features
  heatIn <- dataOrg[, colnames(dataOrg) %in% c("TE", bestSel)]
  heatOut <- scale(heatIn, center = TRUE, scale = TRUE)
  heatOut <- heatOut[, apply(heatOut, 2, function(x) !any(is.na(x)))]
  colnames(heatOut) <- namesDf$originalNames[match(colnames(heatOut), namesDf$newNames)]
  colnames(heatOut)[length(colnames(heatOut))] <- "EffectMeasure"
  
  pdf(paste(nameOut,"heatmap.pdf", sep = "_"), useDingbats = F, width = 20, height = 24)
  par(mar = c(8, 5, 8, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.7, cex.main = 1.3, cex.lab = 1.3)
  col <- colorRampPalette(c("blue", "white", "red"))(256)
  gplots::heatmap.2(heatOut,
                    col = col,
                    breaks = c(seq(-5, 5, length = 257)),
                    margins = c(20, 50),
                    key = TRUE,
                    keysize = 0.5,
                    dendrogram = "both",
                    trace = "none",
                    density.info = "none",
                    labCol = colnames(heatOut),
                    key.par = list(cex = 0.9),
                    # Colv         = NULL,
                    cexCol = 1,
                    cexRow = 0.025,
                    key.xlab = "",
                    lhei = c(3, 25),
                    lwid = c(3, 13),
                    na.rm = TRUE,
                    main = ""
                    # RowSideColors=
  )
  dev.off()
}
