uorf_analysis <- function(ads,
                          regulation,
                          contrast,
                          comparisons=NULL,
                          annot,
                          onlyUTR5=FALSE,
                          startCodon='ATG',#
                          KozakContext='strong',#'adequate1','adequate2','weak','any'
                          selection, #shortest, longest, random (default)
                          geneList=NULL,
                          geneListnames=NULL,
                          geneListcolours=NULL,
                          pdfName=NULL
){
  #Subset annot for only expressed genes
  bg <- row.names(ads@dataP)
  annotBg <- annot[annot$geneID %in% bg,]
  #
  #Select per gene level
  if(selection=='shortest'){
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice(which.min(lenTmp)))
  } else if(selection=='longest'){
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice(which.max(lenTmp)))
  } else {
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice_sample(n = 1))
  }
  #
  if(!isTRUE(onlyUTR5)){
    ext <- as.character(mapply(combSeq, annotBgSel$CDS_seq,annotBgSel$UTR3_seq))
  } else {
    ext <- NULL
  }
  #
  if(KozakContext=='strong'){
    uorfOut <- as.numeric(mapply(calc_uORF, annotBgSel$UTR5_seq, ext=ext, context=paste('[AG][ATGC][ATGC]',startCodon,'G',sep='')))
  } else if(KozakContext=='adequate1'){
    uorfOut <- as.numeric(mapply(calc_uORF, annotBgSel$UTR5_seq, ext=ext, context=paste('[AG][ATGC][ATGC]',startCodon,'[ATC]',sep='')))
  } else if(KozakContext=='adequate2'){
    uorfOut <- as.numeric(mapply(calc_uORF, annotBgSel$UTR5_seq, ext=ext, context=paste('[TC][ATGC][ATGC]',startCodon,'G',sep='')))
  } else if(KozakContext=='weak'){
    uorfOut <- as.numeric(mapply(calc_uORF, annotBgSel$UTR5_seq, ext=ext, context=paste('[TC][ATGC][ATGC]',startCodon,'[ATC]',sep='')))
  } else if(KozakContext=='any'){
    uorfOut <- as.numeric(mapply(calc_uORF, annotBgSel$UTR5_seq, ext=ext, context=paste('[ATGC][ATGC][ATGC]',startCodon,'[ATGC]',sep='')))
  } else {
    stop("Please provide correct Kozak context")
  }
  #
  names(uorfOut) <- annotBgSel$geneID
  #
  #Prepare plotting
  AnotaColours <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],RColorBrewer::brewer.pal(8,"Reds")[c(2,6)],RColorBrewer::brewer.pal(8,"Greens")[c(4,8)], RColorBrewer::brewer.pal(8,"Greens")[c(2,6)],RColorBrewer::brewer.pal(8,"Blues")[c(4,8)])
  names(AnotaColours) <- c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown","bufferingmRNAUp","bufferingmRNADown")
  #
  #Extract all results
  results <- anota2seqGetDirectedRegulations(ads)
  
  res <- vector("list", length = length(regulation))
  for(i in unique(contrast)){
    resTmp <- results[[i]][names(results[[i]]) %in% regulation[which(contrast==i)]]
    res[which(contrast==i)] <- resTmp
  }
  if(!is.null(geneList)){
    res <- append(res, geneList)
  }
  #
  resCol <- list()
  #Load background first
  resCol[[1]] <- as.numeric(uorfOut)
  #
  resProp <- as.numeric()
  #Load background first
  resProp[1] <- length(uorfOut[uorfOut>0])/length(uorfOut)
  #
  for(j in 1:length(res)){
    tmpVal <- as.numeric(uorfOut[names(uorfOut) %in% res[[j]]])
    #
    resProp[[j+1]] <- length(tmpVal[tmpVal>0])/length(tmpVal)
    #
    resCol[[j+1]] <- tmpVal
  }
  #
  if(!is.null(geneList)){
    tmpColour <- append(as.character(AnotaColours[regulation]), geneListcolours)
  } else {
    tmpColour <- as.character(AnotaColours[regulation])
  }
  
  #Plot
  pdf(ifelse(is.null(pdfName),paste('uORFs_',KozakContext,'.pdf',sep=''), paste(pdfName,'_uORFs_',KozakContext,'.pdf',sep='')),width= 8,height=8, useDingbats = F)
  par(mar=c(8,12,5,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
  barplot(resProp,col=c('grey65',tmpColour), xaxt='n', xlab='', ylab='',main="",lwd=1,bty="n",yaxt="n",font=2,ylim=c(0,1+(length(comparisons)*0.1)),space=0)

  axis(side=2, font=2,las=2,lwd=2,at=seq(0,1,0.2),labels = seq(0,1,0.2))
      
  mtext(side=2, line=6, paste('proportion of uORFs in \n',KozakContext,' Kozak context',sep=''), col="black", font=2, cex=1.7,at=0.5)
  #
  if(!is.null(geneList)){
    text(seq(0.5,(length(regulation)+length(geneList)+0.5),1), par("usr")[3] - 0.05, labels=c('background',paste(paste('c',contrast,sep=''),regulation,sep='_'), geneListnames), xpd=NA,cex=0.9,srt=45,adj=1.1)
  } else {
    text(seq(0.5,(length(regulation)+0.5),1), par("usr")[3] - 0.05, labels=c('background',paste(paste('c',contrast,sep=''),regulation,sep='_')), xpd=NA,cex=0.9,srt=45,adj=1.1)
  }
  #
  #Plot stats
  if(!is.null(comparisons)){
    for(j in 1:length(comparisons)){
      compTmp <- comparisons[[j]]+1
      #stats
      pvalTmp <- format(as.numeric(wilcox.test(resCol[[compTmp[1]]], resCol[[compTmp[2]]],alternative='two.sided')[3]),scientific = TRUE,digits=2)
      #
      yposTmp <-  (max(as.numeric(resProp))+0.05)+(j*0.1)
      rect(xleft = compTmp[1]-0.5,xright = compTmp[2]-0.5,ybottom = yposTmp, ytop = yposTmp,lwd=2)
      #
      text((sum(compTmp)/2)-0.5,yposTmp+0.05, pvalTmp ,cex=0.75)
      
    }
  }
  dev.off()
  #
  return(uorfOut)
}

