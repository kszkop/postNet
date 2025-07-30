featureIntegration <- function(ptn,
                               features,
                               lmfeatGroup=NULL,
                               lmfeatGroupColour=NULL,
                               analysis_type,
                               regOnly = TRUE,
                               allFeat = FALSE,
                               useCorel=TRUE,
                               covarFilt = 20,
                               NetModelSel = "omnibus",
                               comparisons = NULL,
                               fdrUni = 0.05,
                               stepP = 0.05,
                               pdfName = NULL) {
  #
  check_ptn(ptn)
  check_features(features)
  if(!is.null(comparisons)){
    if(!check_comparisons(comparisons)){
      stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
    }
    #
    if(length(which(unique(unlist(comparisons))==0))>0 && is.null(ptn_background(ptn))){
      stop(" 0 is always a background, but no background provided")
    }
  }
  check_analysis_type(analysis_type)
  if(!check_logical(regOnly)){
    stop("'regOnly' can only be TRUE or FALSE")
  }
  if(analysis_type=='lm'){
    if(!check_logical(allFeat)){
      stop("'allFeat' can only be TRUE or FALSE")
    }
    if(!check_logical(useCorel)){
      stop("'useCorel' can only be TRUE or FALSE")
    }
    if(!check_number(covarFilt)){
      stop("'covarFilt' can only be a numerical value")
    }
    if(!is_valid_NetModelSel(NetModelSel)){
      stop("'NetModelSel' has to be not null and one one of the: omnibus, adjusted")
    }
  }
  if(analysis_type=='rf'){
    if(is.null(comparisons)){
      stop('Please provide desired comparisons')
    }
  }
  #
  ptn <- prepFeatures(ptn, features)
  
  dataTmp <- ptn_features(ptn)
  
  
  effTmp <- ptn_effect(ptn)
  colnames(dataTmp) <- c(paste("a", seq(1, ncol(dataTmp), 1), sep = ""))
  
  dataTmp$effM <- effTmp[match(row.names(dataTmp),names(effTmp))]
  namesDf <- data.frame(originalNames = colnames(ptn_features(ptn))[1:ncol(dataTmp)-1], newNames = colnames(dataTmp)[1:ncol(dataTmp)-1], stringsAsFactors = F)
  #
  resOut <- resQuant(qvec = ptn_effect(ptn), ptn = ptn)
  #
  #fiOut <- new("postNetFeatureIntegration",
  #             lm = NULL,
  #             rf = NULL,
  #             featureMap = NULL)
  
  ######
  if (analysis_type == 'lm'){
    #
    if(!is.null(lmfeatGroup)){
      check_lmfeatGroup(lmfeatGroup, ncol(dataTmp)-1)
      names(lmfeatGroup) <- colnames(dataTmp)[1:ncol(dataTmp)-1]
      
      if(is.null(lmfeatGroupColour)){
        lmfeatGroupColourOut <- colourAssign(group = lmfeatGroup, colours = lmfeatGroupColour)
      } else {
        check_lmfeatGroupColour(lmfeatGroupColour, lmfeatGroup)
        
        lmfeatGroupColourOut <- colourAssign(group = lmfeatGroup, colours = lmfeatGroupColour)
      }
    }
    #
    if (isTRUE(regOnly)){
      #
      compOut <- list()
      for (i in 1:length(comparisons)) {
        coloursTmp <- ptn_colours(ptn)
        if (names(resOut)[1] == 'background') {
          compTmp <- comparisons[[i]] + 1
          coloursTmp <- c('grey75',coloursTmp)[compTmp]
        } else {
          compTmp <- comparisons[[i]]
          coloursTmp <- coloursTmp[compTmp]
        }
        listSel <- c(names(resOut[[compTmp[1]]]), names(resOut[[compTmp[2]]]))
        dataTmpSel <- dataTmp[row.names(dataTmp) %in% listSel, ]
        
        nameOut <- ifelse(is.null(pdfName), paste('lm', paste(names(resOut)[compTmp], collapse = '_'), sep='_'), paste(pdfName, 'lm', paste(names(resOut)[compTmp], collapse = '_'), sep='_'))
        #
        lmOut <- runLM(dataIn = dataTmpSel, namesDf = namesDf, allFeat = allFeat, useCorel = useCorel, covarFilt=covarFilt, nameOut = nameOut, NetModelSel = NetModelSel, coloursIn=coloursTmp,lmfeatGroup=lmfeatGroup,lmfeatGroupColour=lmfeatGroupColourOut, fdrUni = fdrUni, stepP = stepP)
        compOut[[paste(names(resOut)[compTmp], collapse='_')]] <- lmOut
        
        #fiOut@lm[[paste(names(resOut)[compTmp], collapse='_')]] <- lmOut
        
        bestSel <- names(lmOut@selectedFeatures)
        
        for (feat in bestSel) {
          #
          featTmp <- namesDf[namesDf$originalNames == feat, ]$newNames
          #
          set <- dataTmpSel[,colnames(dataTmpSel) %in% c(featTmp,'effM')]
          #
          set1 <- names(resOut[[compTmp[1]]])
          setSel1 <- set[row.names(set) %in% set1,]
          set2 <- names(resOut[[compTmp[2]]])
          setSel2 <- set[row.names(set) %in% set2,]
          #
          plotScatterInd(set1=setSel1, set2=setSel2, orgName=feat, coloursIn=coloursTmp, nameOut=nameOut)
        }
        #compOut[i] <- paste(names(resOut)[compTmp], collapse='_') 
      }
      #fiOut@lm <- compOut
    } else {
      #fiOut@comparisons <- 'allData'
      #
      coloursTmp <- c('salmon','skyblue')
      lmOut <- runLM(dataIn = dataTmp, namesDf = namesDf, allFeat = allFeat, useCorel = useCorel,  covarFilt=covarFilt, nameOut = pdfName, NetModelSel = NetModelSel, coloursIn=coloursTmp,lmfeatGroup=lmfeatGroup,lmfeatGroupColour=lmfeatGroupColourOut)
      #fiOut@lm[['allData']] <- lmOut
      compOut <- lmOut
      #
      bestSel <- names(lmOut@selectedFeatures)
      
      nameOut <- ifelse(is.null(pdfName), 'lm_allData', paste(pdfName, 'lm_allData', sep='_'))
      
      
      for (feat in bestSel) {
        #
        featTmp <- namesDf[namesDf$originalNames == feat, ]$newNames
        #
        set <- dataTmp[,colnames(dataTmp) %in% c(featTmp,'effM')]
        #
        #set1 <- names(resOut[[compTmp[1]]])
        #setSel1 <- set[row.names(set) %in% set1,]
        #set2 <- names(resOut[[compTmp[2]]])
        #setSel2 <- set[row.names(set) %in% set2,]
        
        plotScatterInd(set1=set, set2=NULL, orgName=feat, coloursIn='grey75', nameOut=nameOut)
      }
    }
    ptn@analysis@featureIntegration[['lm']] <- compOut
  } else if (analysis_type == "rf") {
    dataTmpReg <- dataTmp[, colnames(dataTmp) != "effM"]
    colnames(dataTmpReg) <- namesDf$originalNames[match(colnames(dataTmpReg), namesDf$newNames)]
    
    #
    compOut <- list()
    for (i in 1:length(comparisons)) {
      coloursTmp <- ptn_colours(ptn)
      if (names(resOut)[1] == 'background') {
        compTmp <- comparisons[[i]] + 1
        coloursTmp <- c('grey75',coloursTmp)[compTmp]
      } else {
        compTmp <- comparisons[[i]]
        coloursTmp <- coloursTmp[compTmp]
      }
      #regTmp <- names(resOut)[compTmp]
      nameOut <- ifelse(is.null(pdfName), paste('randomForest', paste(names(resOut)[compTmp], collapse = '_'), sep='_'), paste(pdfName, 'randomForest', paste(names(resOut)[compTmp], collapse = '_'), sep='_'))
      
      dataTmpSel <- dataTmpReg
      dataTmpSel$reg <- NA
      for(j in 1:2) {
        if(length(compTmp) != 2){
          stop('There is something wrong with comparisons')
        }
        cTmp <- names(resOut[[compTmp[j]]])
        regTmp <- c('A','B')
        dataTmpSel$reg[row.names(dataTmpSel) %in% cTmp] <- regTmp[j]
      }
      dataTmpSel <- dataTmpSel[!is.na(dataTmpSel$reg),]
      dataTmpSel$reg <- as.factor(dataTmpSel$reg)
      #
      # Split to Train and Valid ( 70:30)
      train <- sample(nrow(dataTmpSel), 0.7 * nrow(dataTmpSel), replace = FALSE)
      TrainSet <- dataTmpSel[train, ]
      TrainSet$reg <- as.factor(TrainSet$reg)
      ValidSet <- dataTmpSel[-train, ]
      ValidSet$reg <- as.factor(ValidSet$reg)
      # run model on training set
      model1 <- randomForest::randomForest(reg ~ ., data = TrainSet, importance = TRUE, ntree = 500)
      #row.names(model1$importance) <- namesDf$originalNames[match(row.names(model1$importance), namesDf$newNames)]
      #row.names(model1$importanceSD) <- namesDf$originalNames[match(row.names(model1$importanceSD), namesDf$newNames)]
      
      #
      model1Imp <- Boruta::Boruta(reg ~ ., data = TrainSet, doTrace = 0, maxRuns = 500, pValue = 0.001)
      # selecct important once
      featComf <- row.names(Boruta::attStats(model1Imp))[which(as.character(Boruta::attStats(model1Imp)[, 6]) == "Confirmed")]
      #
      pdf(paste(nameOut, "featureImportance.pdf", sep = "_"), width = 8, height = 8, useDingbats = F)
      par(mar = c(10, 5, 3, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
      plot(model1Imp, las = 2, xlab = "", ylab = "", yaxt = "n", xaxt = "n", pch = 20)
      mtext(side = 1, line = 9, "Features", col = "black", font = 2, cex = 1.2)
      mtext(side = 2, line = 3, "Importance (Z-score)", col = "black", font = 2, cex = 1.2)
      axis(side = 2, seq(0, roundNice(max(Boruta::attStats(model1Imp)[, 4]), direction='up'), 10), font = 2, lwd = 2, las = 2, cex = 0.75)
      #
      tmp <- lapply(1:ncol(model1Imp$ImpHistory), function(i) model1Imp$ImpHistory[is.finite(model1Imp$ImpHistory[, i]), i])
      names(tmp) <- colnames(model1Imp$ImpHistory)
      tmpNames <- names(sort(sapply(tmp, median)))
      addNames <- c("shadowMin", "shadowMax", "shadowMean")
      
      #tmpNames <- c(namesDf$originalNames, addNames)[match(tmpNames, c(namesDf$newNames, addNames))]
      coloursN <- rep("black", length(tmpNames))
      coloursN[tmpNames %in% addNames] <- "firebrick1"
      axis(side = 1, at = 1:length(tmpNames), labels = F, font = 2, lwd = 2, las = 2, cex.axis = 0.5,tck=-0.005)
      text(1:length(tmpNames), par("usr")[3] - 1.05, labels = tmpNames, col = coloursN, srt = 45, adj = 1, cex = 0.55, xpd = NA)
      dev.off()
      
      # rerun model with only confirmed ones
      TrainSet <- TrainSet[, colnames(TrainSet) %in% c(featComf, "reg")]
      ValidSet <- ValidSet[, colnames(ValidSet) %in% c(featComf, "reg")]
      
      # run model on training set
      model2 <- randomForest::randomForest(reg ~ ., data = TrainSet, importance = TRUE, ntree = 500)
      #
      varImpIn <- sort(randomForest::importance(model2)[, 3], decreasing = T)
      #names(varImpIn) <- namesDf$originalNames[match(names(varImpIn), namesDf$newNames)]
      #
      pdf(paste(nameOut,"FinalModel.pdf", sep = "_"), width = 16, height = 8, useDingbats = F)
      par(mfrow = c(1, 2), mar = c(9, 5, 10, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
      colDot <- rep("black", length(randomForest::importance(model2)[, 3]))
      #colDot[which(names(sort(randomForest::importance(model2)[, 3], decreasing = F)) %in% featComf)] <- "#B0F2BC"
      dotchart(sort(randomForest::importance(model2)[, 3], decreasing = F), cex = 0.75, col = colDot, labels = names(sort(varImpIn, decreasing = F)), xlab = "", xaxt = "n", frame.plot = FALSE, pch = 16)
        
      axis(side = 1, seq(0, roundNice(max(varImpIn),direction='up'), 5), font = 2, lwd = 2)
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
      
      pdf(paste(nameOut,'pred_rocr.pdf', sep='_'),width=8,height=8, useDingbats = F)
      plot(predOut, main = paste("ROC Curve for Random Forest \n AUC: ", round(auc@y.values[[1]], 3), sep = ""), col = "firebrick1", lwd = 3, xlab = "", ylab = "", )
      abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")
      
      mtext(side = 1, line = 4, "False positive rate", col = "black", font = 2, cex = 1.2)
      mtext(side = 2, line = 3, "True positive rate", col = "black", font = 2, cex = 1.2)
      text(0.8, 0.2, font = 2, cex = 1.7, paste("Sensitivity: ", round(caret::confusionMatrix(predValidc, ValidSet$reg)[[4]][1], 2), sep = ""))
      text(0.8, 0.1, font = 2, cex = 1.7, paste("Specificity: ", round(caret::confusionMatrix(predValidc, ValidSet$reg)[[4]][2], 2), sep = ""))
      dev.off()
      
      rfOut <- new("postNetFeatureIntegration_rf",
                   preModel = model1,
                   borutaModel = model1Imp,
                   finalModel = model2,
                   selectedFeatures = varImpIn)
      #
      compOut[[paste(names(resOut)[compTmp], collapse='_')]] <-  rfOut
      #fiOut@rf[[paste(names(resOut)[compTmp], collapse='_')]] <- rfOut
    
      bestSel <- names(rfOut@selectedFeatures)
        
      for (feat in bestSel) {
        #
        featTmp <- namesDf[namesDf$originalNames == feat, ]$newNames
        #
        set <- ptn_effect(ptn)
        set <- dataTmp[,colnames(dataTmp) %in% c(featTmp,'effM')]
        #
        set1 <- names(resOut[[compTmp[1]]])
        setSel1 <- set[row.names(set) %in% set1,]
        set2 <- names(resOut[[compTmp[2]]])
        setSel2 <- set[row.names(set) %in% set2,]
        #
        plotScatterInd(set1=setSel1, set2=setSel2, orgName=feat, coloursIn=coloursTmp, nameOut=nameOut)
      }
      #compOut[i] <- paste(names(resOut)[compTmp], collapse='_') 
    }
    #fiOut@rf <- compOut
    ptn@analysis@featureIntegration[['rf']] <- compOut
  } else {
    stop("Please provide correct type: lm for linear regression or rf for random forest")
  }
  return(ptn)
}

