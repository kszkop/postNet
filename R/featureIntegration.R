#function that do that
featureIntegration <- function(ads, 
                               contrast, 
                               RegMode='translation', #to choose translatedmRNA, totalmRNA, translation, buffering
                               features, 
                               regOnly=TRUE,
                               pdfName=NULL){
  #
  scOut <- anota2seq::anota2seqGetOutput(ads, output = "singleDf", selContrast = contrast, getRVM = TRUE)
  #
  TE <-  scOut[,grepl(paste(RegMode,'apvEff',sep='.'),colnames(scOut))]
  names(TE) <- scOut$identifier
  
  ##Rename colnames but save original names
  featureName <- names(features)
  #
  features <- append(features, list(TE))
  featPos <- length(features)
  features <- unname(features)
  #
  tmpDf <- data.frame(t(plyr::ldply(features, rbind)))
  tmpDf <- na.omit(tmpDf)
  #
  #TE <- tmpDf[,ncol(tmpDf)]
  #
  dat <- tmpDf#[,1:(ncol(tmpDf))-1]
  colnames(dat) <- c(paste('a',seq(1,length(featureName),1),sep=''),'TE')
  #
  if(isTRUE(regOnly)){
    resTmp <- anota2seqGetDirectedRegulations(ads)[[contrast]]
    #
    if(RegMode=='translation'| RegMode=='translatedmRNA'){  
      listSel <- as.character(unlist(resTmp[c(1,2)]))
    } else if (RegMode=='buffering'){
      listSel <- as.character(unlist(resTmp[c(3,4)]))
    } else if (RegMode=='totalmRNA'){
      listSel <- as.character(unlist(resTmp[c(5,6)]))
    }
    dat <- dat[row.names(dat) %in% listSel,]
  }
  #
  #
  namesDf <- data.frame(originalNames=featureName,newNames=colnames(dat)[-featPos], stringsAsFactors = F)
  #
  #Collect all Fvalues for all models
  fval <- list()
  pval <- list()
  #
  #Start with univariate
  models <- lapply(colnames(dat)[-featPos], function(x) {
    anova(lm(substitute(TE ~ i, list(i = as.name(x))), data = dat))
  })
  #Collect fvalues
  fvalTmp <- sapply(models, function(x) x[1,4])
  names(fvalTmp) <- colnames(dat)[-featPos]
  #if na put 0
  fvalTmp[is.na(fvalTmp)] <- 0
  #
  fval[[1]] <- fvalTmp
  
  #Collect pvalues
  pvalTmp <- sapply(models, function(x) x[1,5])
  names(pvalTmp) <- colnames(dat)[-featPos]
  #if na put 0
  pvalTmp[is.na(pvalTmp)] <- 1
  #
  pval[[1]] <- pvalTmp
  
  #Collect features
  bestTmp <- which.max(sapply(models, function(x) x[1,4]))
  outTmp <- which(sapply(models, function(x) x[1,5])>0.05 | is.na(sapply(models, function(x) x[1,5])))
  #
  bestSel <- colnames(dat)[-featPos][bestTmp]
  outSel <- colnames(dat)[-featPos][outTmp]
  
  #
  #cat('\n','Highest Fval: ', namesDf$originalNames[namesDf$newNames %in% tmpIn[bestTmp]])
  #if(length(outTmp) > 0){
  #  cat('    No longer significant: ', namesDf$originalNames[namesDf$newNames %in% tmpIn[outTmp]], '\n')
  #} else {
  #  cat('    Non removed ','\n')
  #}
  #
  i <- 1
  #Start with multivariate
  while(length(colnames(dat)[-featPos][!colnames(dat)[-featPos] %in% c(bestSel,outSel)])>0){
    #Counter
    i <- i + 1
    #
    tmpIn <- colnames(dat)[-featPos][!colnames(dat)[-featPos] %in% c(bestSel,outSel)]
    models2 <- list()
    #
    x_sel <- paste(bestSel, collapse=" + ")
    
    for(j in 1:(length(colnames(dat)[-featPos])-length(c(bestSel,outSel)))){
      varx <-  colnames(dat)[-featPos][!colnames(dat)[-featPos] %in% c(bestSel,outSel)][j]
      design <- as.formula(paste(paste('TE', paste(x_sel, collapse=" + "), sep="~"),varx,sep ='+'))
      models2[[j]] <-  anova(lm(design, data = dat))
    }
    #Collect fvalues
    fvalTmp <- sapply(models2, function(x) x[nrow(x)-1,4])
    names(fvalTmp) <- colnames(dat)[-featPos][!colnames(dat)[-featPos] %in% c(bestSel,outSel)]
    
    fval[[i]] <- fvalTmp
    
    #Collect pvalues
    pvalTmp <- sapply(models2, function(x) x[nrow(x)-1,5])
    names(pvalTmp) <- colnames(dat)[-featPos][!colnames(dat)[-featPos] %in% c(bestSel,outSel)]
    
    pval[[i]] <- pvalTmp
    #
    bestTmp <- which.max(sapply(models2, function(x) x[nrow(x)-1,4]))
    outTmp <- which(sapply(models2, function(x) x[nrow(x)-1,5])>0.05)
    #
    if(length(outTmp)>0){
      if(!bestTmp %in% outTmp){
        bestSel <- append(bestSel,tmpIn[bestTmp])
      }
    } else {
      bestSel <- append(bestSel,tmpIn[bestTmp])
    }
    outSel <- append(outSel,tmpIn[outTmp])
    #
    #cat('\n','Highest Fval: ', namesDf$originalNames[namesDf$newNames %in% tmpIn[bestTmp]])
    #if(length(outTmp) > 0){
    #  cat('    No longer significant: ', namesDf$originalNames[namesDf$newNames %in% tmpIn[outTmp]], '\n')
    #} else {
    #  cat('    Non removed ','\n')
    #}
  }
  #format
  tmp <- t(plyr::ldply(fval, rbind))
  tmp <- tmp[c(bestSel, rev(outSel)),]
  #
  tmp_pval <- t(plyr::ldply(pval, rbind))
  tmp_pval <- tmp_pval[c(bestSel, rev(outSel)),]
  #
  tb1Out <- round(apply(tmp, 2 ,as.numeric),2)
  row.names(tb1Out) <- namesDf$originalNames[match(row.names(tmp), namesDf$newNames)]
  tb1Out[is.na(tb1Out)] <- ""
  
  #PLot table
  colours <- matrix("white", nrow(tb1Out), ncol(tb1Out))
  for(k in 1:ncol(colours)){
    colours[k,k] <- 'seagreen1'
  }
  colours[which(tmp_pval>0.05)] <- 'red3'
  
  tt1 <- gridExtra::ttheme_default(core=list(fg_params=list(fontface=c(rep("plain",ncol(tb1Out)))), bg_params = list(fill =colours , col="black")))
  tg1 <- gridExtra::tableGrob(tb1Out, theme = tt1)
  
  ##Calculate % varaince explained and plot output
  #Final model
  varExpl <- anova(lm(as.formula(paste('TE', paste(bestSel, collapse=" + "), sep="~")), data = dat))
  varExplss <- varExpl$"Sum Sq"
  #Outtable
  tb2out <- data.frame(Features=namesDf$originalNames[match(bestSel, namesDf$newNames)],Pvalue=format(varExpl$`Pr(>F)`[1:length(bestSel)],scientific = T,digits = 2),VarianceExplained=round((varExplss/sum(varExplss)*100)[1:length(bestSel)],2))
  #
  tg2 <- gridExtra::tableGrob(tb2out,rows = NULL)
  
  pdf(ifelse(is.null(pdfName),paste(RegMode,'featureIntegration.pdf',sep='_'), paste(pdfName, RegMode,'featureIntegration.pdf',sep='_')),width= 24,height=6+(length(features)-1)/2, useDingbats = F)
  gridExtra::grid.arrange(tg1, tg2, ncol = 2, nrow = 1)
  dev.off()
}
