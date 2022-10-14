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
  step1expl <- step1expl[-presel]
  step1pval <- step1pval[-presel]
  step1pval_fdr <- step1pval_fdr[-presel]
  #
  #remove this features
  if(!isTRUE(allFeat)){
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
  #colour warnings
  for(k in 2:nrow(colours)){
    trow <- as.numeric(na.omit(as.numeric(tb1Out[k,])))
    #check whihch one is more than 50% drop
    tsel <- which((lag(trow)/trow)>2)
    #
    colours[k,tsel] <- '#FDE0C5'
  }
  #
  for(k in 1:ncol(colours)){
    colours[k,k] <- '#B0F2BC'
  }
  #
  colours[which(tmp_pval>0.05)] <- '#B14D8E'
  
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
  step1sel <- namesDf$newNames[-presel]
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
  
  #
  #Prepare plotting
  AnotaColours <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],RColorBrewer::brewer.pal(8,"Reds")[c(2,6)],RColorBrewer::brewer.pal(8,"Greens")[c(4,8)], RColorBrewer::brewer.pal(8,"Greens")[c(2,6)],RColorBrewer::brewer.pal(8,"Blues")[c(4,8)])
  names(AnotaColours) <- c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown","bufferingmRNAUp","bufferingmRNADown")
  #
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
        set1 <- dat[row.names(dat) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[1]]), colnames(dat)==feat]
        set2 <- dat[row.names(dat) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[2]]), colnames(dat)==feat]
        set_te1 <- dat$TE[row.names(dat) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[1]])]
        set_te2 <- dat$TE[row.names(dat) %in% as.character(resTmp[grepl(regulation,names(resTmp))][[2]])]
        col1 <- as.character(AnotaColours[grepl(regulation,names(AnotaColours))][1])
        col2 <- as.character(AnotaColours[grepl(regulation,names(AnotaColours))][2])
      } else {
        set1 <- dat[row.names(dat) %in% as.character(unlist(resTmp[c(1)])), colnames(dat)==feat]
        set2 <- dat[row.names(dat) %in% as.character(unlist(resTmp[c(2)])), colnames(dat)==feat]
        set_te1 <- dat$TE[row.names(dat) %in% as.character(unlist(resTmp[c(1)]))]
        set_te2 <- dat$TE[row.names(dat) %in% as.character(unlist(resTmp[c(2)]))]
        col1 <- geneListcolours[1]
        col2 <- geneListcolours[2] 
      }
      #
    } else {
      #
      set <- dat[,colnames(dat)==feat]
      set_te <- dat$TE
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
}
