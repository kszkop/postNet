featureIntegration <- function(ads,
                               contrast,
                               RegMode='translation', #to choose translatedmRNA, totalmRNA, translation, buffering
                               geneList=NULL,
                               geneListnames=NULL,
                               geneListcolours=NULL,
                               features,
                               regOnly=TRUE,
                               regulation=NULL, #required if regOnly is TRUE, to choose as it is in anota2seq
                               allFeat=TRUE,
                               pdfName=NULL,
                               type #lm (linear model), rf (random forest)
){
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
  #TE <- tmpDf[,ncol(tmpDf)]
  #
  dat <- tmpDf#[,1:(ncol(tmpDf))-1]
  colnames(dat) <- c(paste('a',seq(1,length(featureName),1),sep=''),'TE')
  #save original names and data
  namesDf <- data.frame(originalNames=featureName,newNames=colnames(dat)[-length(colnames(dat))], stringsAsFactors = F)
  dataOrg <- dat
  #
  if(is.null(geneList)){
    #
    if(isTRUE(regOnly)){
      resTmp <- anota2seqGetDirectedRegulations(ads)[[contrast]]
      #
      listSel <- as.character(unlist(resTmp[grepl(regulation,names(resTmp))]))
      #
      dat <- dat[row.names(dat) %in% listSel,]
    }
  } else {
    #
    resTmp <- geneList
    names(resTmp) <- geneListnames
    #
    listSel <- as.character(unlist(geneList))
    dat <- dat[row.names(dat) %in% listSel,]
  }
  if(type=='rf'){
    if(isTRUE(regOnly)){
      #Change numeric TE to categorical 
      dat$reg <- ifelse(dat$TE>0, 2, 1)
    } else {
      if(is.null(geneList)){
        dat$reg <- 0
        dat$reg[row.names(dat) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[1]])] <- 2
        dat$reg[row.names(dat) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[2]])] <- 1
      } else {
        dat$reg <- 0
        dat$reg[row.names(dat) %in% as.character(unlist(resTmp[c(1)]))] <- 2
        dat$reg[row.names(dat) %in% as.character(unlist(resTmp[c(2)]))] <- 1
      }
    }
    dat <- dat[,colnames(dat) !='TE']
    #Split to Train and Valid ( 70:30)
    train <- sample(nrow(dat), 0.7*nrow(dat), replace = FALSE)
    TrainSet <- dat[train,]
    TrainSet$reg <- as.factor(TrainSet$reg)
    ValidSet <- dat[-train,]
    ValidSet$reg <- as.factor(ValidSet$reg)
    #run model on training set
    model1 <- randomForest::randomForest(reg ~ ., data = TrainSet, importance = TRUE, ntree = 500)
    
    #Plot importance
    #also apply to find relevant features
    model1Imp <- Boruta::Boruta(reg ~ ., data = TrainSet, doTrace = 0, maxRuns = 500,pValue = 0.001)
    #selecct important once
    featComf <- row.names(Boruta::attStats(model1Imp))[which(as.character(Boruta::attStats(model1Imp)[,6]) ==  "Confirmed")]
    #featComf <- namesDf$originalNames[match(featComf, namesDf$newNames)]
    #pdf(ifelse(is.null(pdfName),paste(RegMode,'featureSel_randomForest.pdf',sep='_'), paste(pdfName, RegMode,'featureSel_randomForest.pdf',sep='_')),width=8,height=8, useDingbats = F)
    #plot(model1Imp, las = 2, cex.axis = 0.7)
    #dev.off()
    
    
    varImpIn <- sort(randomForest::importance(model1)[,3],decreasing=T)
    #
    pdf(ifelse(is.null(pdfName),paste(RegMode,'randomForest.pdf',sep='_'), paste(pdfName, RegMode,'randomForest.pdf',sep='_')),width=16,height=8, useDingbats = F)
    par(mfrow=c(1,2),mar=c(9,5,10,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.3, cex.main=1.7,cex.lab=1)
    colDot <- rep('black',length(randomForest::importance(model1)[,3]))
    colDot[which(names(sort(randomForest::importance(model1)[,3],decreasing=F)) %in% featComf)] <- '#B0F2BC'
    dotchart(sort(randomForest::importance(model1)[,3],decreasing=F),cex=0.75, col=colDot,labels = namesDf$originalNames[match(names(sort(randomForest::importance(model1)[,3],decreasing=F)), namesDf$newNames)],xlab='',xaxt='n',frame.plot=FALSE,pch=16)
    
    axis(side=1,seq(0,roundUpNice(max(varImpIn)),5), font=2,lwd=2)
    mtext(side=1, line=4, "Feature Importance \n (Mean Decrease Accuracy)", col='black', font=2,cex=1.2)
    
    conf <- model1$confusion[,-ncol(model1$confusion)]
    oob <- (1 - (sum(diag(conf))/sum(conf)))*100
    mtext(side=3, line=3, paste('OOB estimate of  error rate: ',round(oob,2),sep=''), col="black", font=2, cex=1.7)
    
    #Mean Decrease Accuracy - How much the model accuracy decreases if we drop that variable.
    #Mean Decrease Gini - Measure of variable importance based on the Gini impurity index used for the calculation of splits in trees.
    predValidc <- stats::predict(model1, ValidSet, type = "class")
    
    #Accuracy
    #as.numeric(confusionMatrix(predValidc , ValidSet$reg)[[3]][1])
    #Sensitivity
    #confusionMatrix(predValidc , ValidSet$reg)[[4]][1]
    #Specificity
    #confusionMatrix(predValidc , ValidSet$reg)[[4]][2]
    #
    predValid <- stats::predict(model1, ValidSet, type = "prob")
    #
    perf = ROCR::prediction(predValid[,2], as.numeric(ValidSet$reg))
    # 1. Area under curve
    auc = ROCR::performance(perf, "auc")
    # 2. True Positive and Negative Rate
    predOut = ROCR::performance(perf, "tpr","fpr")
    # 3. Plot the ROC curve
    plot(predOut,main=paste("ROC Curve for Random Forest \n Accuracy: ",round(as.numeric(caret::confusionMatrix(predValidc , ValidSet$reg)[[3]][1]),3),sep=''),col='firebrick1',lwd=3,xlab='',ylab='',)
    abline(a=0,b=1,lwd=2,lty=2,col="gray")
    
    mtext(side=1, line=4, 'False positive rate', col='black', font=2,cex=1.2)
    mtext(side=2, line=3, 'True positive rate', col="black", font=2, cex=1.2)
    text(0.8,0.2, font=2,cex=1.7,paste('Sensitivity: ',round(caret::confusionMatrix(predValidc , ValidSet$reg)[[4]][1],2),sep=''))
    text(0.8,0.1, font=2,cex=1.7,paste('Specificity: ',round(caret::confusionMatrix(predValidc , ValidSet$reg)[[4]][2],2),sep=''))
    
    #axis(side=1,seq(0,1,0.25), font=2,lwd=2)
    #axis(side=2,seq(0,1,0.25), font=2,las=2,lwd=2)
    dev.off()
    
  } else if (type=='lm'){
    
    #Pre-run to remove features that are not significnt for univariate models
    #
    #Start with univariate
    models_prerun <- lapply(colnames(dat)[-featPos], function(x) {
      anova(lm(substitute(TE ~ i, list(i = as.name(x))), data = dat))
    })
    #Select the ones not significant after first step
    presel <- which(sapply(models_prerun , function(x) x[1,5])>0.05 | is.na(sapply(models_prerun , function(x) x[1,5])))
    #
    step1expl <- round(sapply(models_prerun , function(x) (x[1,2]/(sum(x[1,2], x[2,2])))*100),2)
    names(step1expl) <- featureName
    #
    step1pval <- sapply(models_prerun , function(x) (x[1,5]))
    names(step1pval) <- featureName
    #add fdr test
    step1pval_fdr <- p.adjust(step1pval)
    names(step1pval_fdr) <- featureName
    
    #remove not signigicant
    if(length(presel)>0){
      step1expl <- step1expl[-presel]
      step1pval <- step1pval[-presel]
      step1pval_fdr <- step1pval_fdr[-presel]
    }
    #
    #remove this features
    if(!isTRUE(allFeat) & length(presel)>0){
      dat <- dat[,-presel]
      featureName <- featureName[-presel]
    }
    featPos <- length(colnames(dat))
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
    #
    colnames(tb1Out) <- paste('step',seq(1,ncol(tb1Out),1),sep='')
    #
    #PLot table
    colours <- matrix("white", nrow(tb1Out), ncol(tb1Out))
    linkIn <- matrix(NA, nrow(tb1Out), ncol(tb1Out))
    #colour warnings
    for(k in 2:nrow(colours)){
      trow <- as.numeric(na.omit(as.numeric(tb1Out[k,])))
      #check whihch one is more than 50% drop
      tsel <- which((lag(trow)/trow)>2)
      tperc <- trow/lag(trow)*100
      #
      colours[k,tsel] <- '#FDE0C5'
      linkIn[k,1:length(tperc)] <- tperc 
    }
    #
    for(k in 1:ncol(colours)){
      colours[k,k] <- '#B0F2BC'
    }
    #
    colours[which(tmp_pval>0.05)] <- '#B14D8E'
    linkIn[which(linkIn>50)] <- NA
    row.names(linkIn) <- row.names(tmp)
    #
    tt1 <- gridExtra::ttheme_default(core=list(fg_params=list(fontface=c(rep("plain",ncol(tb1Out)))), bg_params = list(fill =colours , col="black")))
    tg1 <- gridExtra::tableGrob(tb1Out, theme = tt1)
    
    ##Calculate % varaince explained and plot output
    #Final model
    varExpl <- anova(lm(as.formula(paste('TE', paste(bestSel, collapse=" + "), sep="~")), data = dat))
    varExpldepend <- numeric()
    for(i in 1:length(bestSel)){
      tmpFeatI <- bestSel[i]
      #variance expl
      varExpldepend[i] <- round((varExpl[i,2]/sum(varExpl$`Sum Sq`))*100,2)
    }
    names(varExpldepend) <- bestSel
    ##CAalculate dependent variance explained loop through the sel features moving each to the end
    #only features in final model
    varExplIndepend <- numeric()
    for(i in 1:length(bestSel)){
      tmpFeat <- bestSel[i]
      restFeat <- bestSel[-i]
      #model
      design <- as.formula(paste(paste('TE', paste(restFeat, collapse=" + "), sep="~"),tmpFeat,sep ='+'))
      tmpM <-  anova(lm(design, data = dat))
      
      #variance expl
      varExplIndepend[i] <- round((tmpM[nrow(tmpM)-1,2]/sum(tmpM$`Sum Sq`))*100,2)
    }
    names(varExplIndepend) <- bestSel
    
    #features significant after step1
    varExplIndepend2 <- numeric()
    #features aftet 1st
    if(length(presel)>0){
      step1sel <- namesDf$newNames[-presel]
    } else {
      step1sel <- namesDf$newNames
    }
    for(i in 1:length(step1sel)){
      tmpFeat2 <- step1sel[i]
      restFeat2<- step1sel[-i]
      #model
      design2 <- as.formula(paste(paste('TE', paste(restFeat2, collapse=" + "), sep="~"),tmpFeat2,sep ='+'))
      tmpM2 <-  anova(lm(design2, data = dat))
      
      #variance expl
      varExplIndepend2[i] <- round((tmpM2[nrow(tmpM2)-1,2]/sum(tmpM2$`Sum Sq`))*100,2)
    }
    names(varExplIndepend2) <- step1sel
    
    #Outtable
    tb2out <- data.frame(Features=namesDf$originalNames[match(bestSel, namesDf$newNames)],Pvalue=format(varExpl$`Pr(>F)`[1:length(bestSel)],scientific = T,digits = 2),VarianceExplained_Omnibus=as.numeric(varExpldepend), VarianceExplained_Adjusted =as.numeric(varExplIndepend))
    
    tg2 <- gridExtra::tableGrob(tb2out,rows = NULL)
    
    tb3out <- data.frame(Features=names(step1expl),Pvalue_Univariate=format(as.numeric(step1pval),scientific = T,digits = 2), FDRvalue_Univariate=format(as.numeric(step1pval_fdr),scientific = T,digits = 2),VarianceExplained_Univariate=as.numeric(step1expl))
    tb3out <- tb3out[with(tb3out,order(-tb3out$VarianceExplained_Univariate)),]
    #
    tg3 <- gridExtra::tableGrob(tb3out,rows = NULL)
    
    pdf(ifelse(is.null(pdfName),paste(RegMode,'featureIntegration.pdf',sep='_'), paste(pdfName, RegMode,'featureIntegration.pdf',sep='_')),width= dim(tg3)[2] + dim(tg1)[2] + dim(tg2)[2]+14.5,height=dim(tg1)[1]/2, useDingbats = F)
    gridExtra::grid.arrange(tg3, tg1, tg2, ncol = 3, nrow = 1,padding=0, top=0, left=0)
    grid::grid.text(paste('Total variance explained: ', sum(as.numeric(varExpldepend)),'%', sep=''), x = grid::unit(0.75, "npc"),  y = grid::unit(0.90, "npc"),gp = grid::gpar(fontsize=15))
    dev.off()
    
    #Outtable
    tb4out <- data.frame(Features=namesDf$originalNames[match(names(varExplIndepend2), namesDf$newNames)],VarianceExplained_IndependentAll=as.numeric(varExplIndepend2))
    tb4out <- tb4out[with(tb4out,order(-tb4out$VarianceExplained_IndependentAll)),]
    #
    tg4 <- gridExtra::tableGrob(tb4out,rows = NULL)
    
    pdf(ifelse(is.null(pdfName),paste(RegMode,'featureIntegration_varexpl_independAll.pdf',sep='_'), paste(pdfName, RegMode,'featureIntegration_varexpl_independAll.pdf',sep='_')),width= dim(tg4)[2] + dim(tg4)[2]+4,height=dim(tg4)[1]/2, useDingbats = F)
    gridExtra::grid.arrange(tg4, ncol = 1, nrow = 1)
    dev.off()
    
    #Plot network
    linkOut <- list()
    for(i in 2:ncol(linkIn)){
      tmpIn <- linkIn[,i]
      #
      tmpOut <- as.numeric(tmpIn)[which(tmpIn>0)]
      if(length(tmpOut)>0){
        names(tmpOut) <- paste(row.names(linkIn)[i-1], names(tmpIn)[which(tmpIn>0)],sep='_')
      }
      #
      linkOut[[i-1]] <- tmpOut
    }
    linkOut <- as.data.frame(unlist(linkOut))
    linkOut <- with(linkOut,cbind(linkOut,reshape2::colsplit(row.names(linkOut),pattern="\\_",names = c('from','to'))))
    rownames(linkOut) <- NULL
    colnames(linkOut)[1] <- 'weight'
    linkOut <- linkOut[,c(2,3,1)] 
    
    ###tb3out
    nodeOut <- tb2out[,c(1,3)]
    nodeOut$ID <- namesDf$newNames[match(nodeOut$Features,namesDf$originalNames)]
    nodeOut <- nodeOut[,c(3,2,1)]
    nodeOut$omnibus <- 1
    
    #other tgat are connected but not significant
    addAll <- data.frame(ID=unique(c(linkOut$from,linkOut$to)), VarianceExplained_Omnibus =0)
    addAll$Features <- namesDf$originalNames[match(addAll$ID,namesDf$newNames)]
    
    addAll <- addAll[!addAll$ID %in% nodeOut$ID,]
    addAll$omnibus <- 2
    #
    nodeOutAll <- rbind(nodeOut,addAll)
    #create igraph object
    net <- igraph::graph.data.frame(linkOut, nodeOutAll, directed=T) 
    #rescale to size ans other attributes
    lsize <-  rescale(igraph::V(net)$VarianceExplained_Omnibus, 0, 100, 0, 50)
    lsize[which(lsize>0)] <- lsize[which(lsize>0)] + 4
    igraph::V(net)$size <- lsize
    lcol <- rep("black",nrow(nodeOutAll))
    lcol[which(nodeOutAll$omnibus==2)] <- "#B14D8E"
    igraph::V(net)$label.color <-  lcol
    igraph::V(net)$label <- wrapNames(igraph::V(net)$Features,8)
    
    #colour nodes as in table
    colrs <- c('#B0F2BC','white')#"#B14D8E")
    igraph::V(net)$color <- colrs[igraph::V(net)$omnibus]
    
    # Set edges width based on weight:
    igraph::E(net)$width <- rescale(igraph::E(net)$weight, 0, 50, 0, 5)
    
    #change arrow size and edge color:
    igraph::E(net)$arrow.size <- .0
    igraph::E(net)$edge.color <- "gray75"
    
    # We can also override the attributes explicitly in the plot:
    pdf(ifelse(is.null(pdfName),paste(RegMode,'network.pdf',sep='_'), paste(pdfName, RegMode,'network.pdf',sep='_')),height=8,width=8, useDingbats = F)
    par(bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.9,cex.main=1.9,cex.lab=1.5)
    #m <- layout(matrix(c(seq(1,nSample+3),rep(c(seq(1,nSample+2),nSample+4),3)), ncol=1), heights=c(1,3))
    plot(net,shape ="sphere",vertex.label.font=2,vertex.label.cex=1,vertex.frame.color="white",layout=layoutCalc(net, n=2))
    
    legend("bottomleft", fill=c('#B0F2BC'), "In Omnibus model",cex=1.3,bty='n',xpd=T,inset=-0.1)
    legend(-1.5,1.5,lwd = c(7,5,3), col="gray75", title=c("Interaction strength"),c("","",""),bty='n',xpd=T)
    legend(-0.8,1.5,pt.cex = c(4,3,2), pch=20,col="gray75", title=c("Variance explained"),c("","",""),bty='n',xpd=T)
    
    dev.off()
  } else {
    stop("Please provide correct type: lm for linear regression or rf for random forest")
  }
  #
  #Prepare plotting
  AnotaColours <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],RColorBrewer::brewer.pal(8,"Reds")[c(2,6)],RColorBrewer::brewer.pal(8,"Greens")[c(4,8)], RColorBrewer::brewer.pal(8,"Greens")[c(2,6)],RColorBrewer::brewer.pal(8,"Blues")[c(4,8)])
  names(AnotaColours) <- c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown","bufferingmRNAUp","bufferingmRNADown")
  #
  if(type=='rf'){
    bestSel <- featComf
  }
  for(feat in bestSel){
    #
    featTmp <- namesDf[namesDf$newNames==feat,]$originalNames
    featTmp <- gsub(' ','_',featTmp)
    #
    if(isTRUE(regOnly)){
      #
      set <- dat[row.names(dat) %in% listSel, colnames(dat)==feat]
      set_te <- dat$TE[row.names(dat) %in% listSel]
      #
      if(is.null(geneList)){
        set1 <- dataOrg[row.names(dataOrg) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[1]]), colnames(dataOrg)==feat]
        set2 <- dataOrg[row.names(dataOrg) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[2]]), colnames(dataOrg)==feat]
        set_te1 <- dataOrg$TE[row.names(dataOrg) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[1]])]
        set_te2 <- dataOrg$TE[row.names(dataOrg) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[2]])]
        col1 <- as.character(AnotaColours[grepl(regulation,names(AnotaColours))][1])
        col2 <- as.character(AnotaColours[grepl(regulation,names(AnotaColours))][2])
      } else {
        set1 <- dataOrg[row.names(dataOrg) %in% as.character(unlist(resTmp[c(1)])), colnames(dataOrg)==feat]
        set2 <- dataOrg[row.names(dataOrg) %in% as.character(unlist(resTmp[c(2)])), colnames(dataOrg)==feat]
        set_te1 <- dataOrg$TE[row.names(dataOrg) %in% as.character(unlist(resTmp[c(1)]))]
        set_te2 <- dataOrg$TE[row.names(dataOrg) %in% as.character(unlist(resTmp[c(2)]))]
        col1 <- geneListcolours[1]
        col2 <- geneListcolours[2] 
      }
      #
    } else {
      #
      set <- dataOrg[,colnames(dataOrg)==feat]
      set_te <- dataOrg$TE
      #
    }
    pdf(ifelse(is.null(pdfName),paste(RegMode,featTmp,'individually.pdf',sep='_'), paste(pdfName, RegMode,featTmp,'individually.pdf',sep='_')),width=8,height=8, useDingbats = F)
    par(mar=c(9,5,5,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.7, cex.main=1.7,cex.lab=1.3)
    
    #
    xlim_max <- ifelse(max(set)>=0, roundUpNice(abs(max(set))),  -roundUpNice(abs(max(set))))
    if(xlim_max < max(set)){
      xlim_max <- ceiling(max(set))
    }
    xlim_min <- ifelse(min(set)>=0, roundUpNice(abs(min(set))),  -roundUpNice(abs(min(set))))
    if(xlim_min > min(set)){
      xlim_min <- floor(min(set))
    }
    ylim_max <- roundUpNice(abs(max(set_te)))
    #
    plot(set,set_te,col='#8A8683',pch=16,cex=1,ylim=c(-ylim_max,ylim_max),xlab='',ylab='',lwd=1,bty="n",xaxt="n",yaxt="n",font=2, xlim=c(xlim_min, xlim_max))
    
    if(isTRUE(regOnly)){
      points(set1, set_te1, pch=16,col=col1)
      points(set2, set_te2, pch=16,col=col2)
    }
    
    #
    mtext(side=2, line=3, RegMode, col="black", font=2, cex=1.7)
    axis(side=2,seq(-ylim_max,ylim_max,2), font=2,las=2,lwd=2)
    
    mtext(side=1, line=4, featTmp, col="black", font=2, cex=1.7,at=(xlim_min+xlim_max)/2)
    axis(side=1,seq(xlim_min,xlim_max,ifelse((xlim_max - xlim_min)/5 >=0 , roundUpNice((xlim_max - xlim_min)/5), -roundUpNice(abs((xlim_max - xlim_min)/5)))), font=2,lwd=2)
    #
    if(length(unique(set))>2 & IQR(set) > 0 & length(unique(set))>3){
      f1 <-predict(smooth.spline(set_te~set))
      lines(f1$x[which(f1$x>xlim_min & f1$x<xlim_max)],f1$y[which(f1$x>xlim_min & f1$x<xlim_max)], col='#AFBADC',lwd=4,lend=2)
      lines(f1$x[which(f1$x>xlim_min & f1$x<xlim_max)],f1$y[which(f1$x>xlim_min & f1$x<xlim_max)], col='black',lwd=1,lend=2,lty=3)
    }
    #
    if(!is.na(as.numeric(coefficients(lm(set_te~set))[2]))){
      plotrix::ablineclip(lm(set_te~set), col='#AFBADC',lwd=4,x1=xlim_min,x2=xlim_max)
      plotrix::ablineclip(lm(set_te~set), col='black',lwd=1,x1=xlim_min,x2=xlim_max)
      
      text((xlim_min+xlim_max)/2,ylim_max, paste('pvalue ',format(as.numeric(cor.test(set,set_te)[3]),scientific=T,digits=3),', r=',round(as.numeric(cor.test(set,set_te)[4]),3),sep=''), bty='n',col='black',cex=1.25,font=2)
    }
    dev.off()
  }
  #heatmap select only important/omnibus features
  heatIn <- dataOrg[,colnames(dataOrg) %in% c('TE',bestSel)]
  heatOut <- scale(heatIn, center = TRUE, scale = TRUE)
  heatOut  <- heatOut[ ,apply(heatOut , 2, function(x) !any(is.na(x)))]
  colnames(heatOut) <-  namesDf$originalNames[match(colnames(heatOut),namesDf$newNames)]
  colnames(heatOut)[length(colnames(heatOut))] <- 'EffectMeasure'
  
  pdf(ifelse(is.null(pdfName),paste(RegMode,featTmp,'heatmap.pdf',sep='_'), paste(pdfName, RegMode,featTmp,'heatmap.pdf',sep='_')), useDingbats = F,width= 20,height=24)
  par(mar=c(8,5,8,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.7, cex.main=1.3,cex.lab=1.3)
  col <- colorRampPalette(c("blue","white", "red"))(256)
  gplots::heatmap.2(heatOut,
                    col          = col,
                    breaks       = c(seq(-5,5,length=257)),
                    margins      = c(20, 50),
                    key =TRUE,
                    keysize      = 0.5,
                    dendrogram   = "row", 
                    trace        = "none", 
                    density.info = "none",
                    labCol       = colnames(heatOut),
                    key.par = list(cex=0.9),
                    Colv         = NULL,
                    cexCol       = 1,
                    cexRow       = 0.025,
                    key.xlab     = "",
                    lhei=c(3,25), 
                    lwid=c(3,13),
                    na.rm=TRUE,
                    main=''
                    #RowSideColors=
  )
  dev.off()
}

