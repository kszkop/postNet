#Function based on given AA/codons calculate freq across all genes
codonCalc  <- function(codonIn, #output of codonUsage function
                       analysis, #whether 'AA' or 'codon'
                       featsel, #vector of selected codons
                       ads,
                       regulation,
                       contrast,
                       comparisons=NULL,
                       unit='count', #option 'freq'
                       pdfName=NULL
){
  
  if(analysis=='codon'){
    codonTmp <- codonIn[codonIn$codon %in% featsel,]
  } else if(analysis=='AA'){
    codonTmp <- codonIn[codonIn$AA %in% featsel,]
  }
  if(unit=='count'){
    tmp <- codonTmp %>% group_by(geneID) %>% summarise(count=sum(codonCount))
    codonCalcOut <- tmp$count
  } else if(unit=='freq'){
    tmp <- codonTmp %>% group_by(geneID) %>% summarise(freq=sum(codonFreq))
    codonCalcOut <- tmp$freq
  }
  names(codonCalcOut) <- tmp$geneID
  
  #
  #Extract all results
  results <- anota2seqGetDirectedRegulations(ads)
  #
  AnotaColours <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],RColorBrewer::brewer.pal(8,"Greens")[c(4,8)], RColorBrewer::brewer.pal(8,"Blues")[c(4,8)])
  names(AnotaColours) <- c("translationUp","translationDown", "mRNAAbundanceUp","mRNAAbundanceDown","bufferingmRNAUp","bufferingmRNADown")
  #subset for desired regulations depending on number of contrasts
  #if(length(unique(contrast))==1){
  #  res <- results[[1]][names(results[[1]]) %in% regulation]
  #} else {
  res <- vector("list", length = length(regulation))
  for(i in unique(contrast)){
    resTmp <- results[[i]][names(results[[i]]) %in% regulation[which(contrast==i)]]
    res[which(contrast==i)] <- resTmp
  }

  #Prepare comparisons per regulation
  for(j in 1:length(comparisons)){
    compTmp <- comparisons[[j]]
    #
    resComp <- res[compTmp]
    
    #Plot
    pdf(ifelse(is.null(pdfName),paste(paste(regulation[compTmp],collapse='_'),'codonUsage.pdf',sep='_'), paste(pdfName,paste(regulation[compTmp],collapse='_'),'codonUsage.pdf',sep='_')),width= 8,height=8, useDingbats = F)
    par(mar=c(5,5,8,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
    plot(ecdf(as.numeric(codonCalcOut)),col='grey45',main='',xlab='',ylab='',verticals=TRUE, do.p=FALSE,lwd=3,bty="n",yaxt="none",font=2, xlim=c(0,as.numeric(quantile(as.numeric(codonCalcOut),0.99))),xaxt="none")
  
    mtext(side=1, line=4, paste('codon usage \n',unit,sep=''), col='black', font=2,cex=1.2)
    mtext(side=2, line=3, 'Fn(x)', col="black", font=2, cex=1.2)
  
    axis(side=1,seq(0,as.numeric(quantile(as.numeric(codonCalcOut),0.99)),1), font=2,lwd=2)
    axis(side=2,seq(0,1,0.2), font=2,las=2,lwd=2)
  
    #
    tableOut <- matrix(NA, nrow= length(regulation[compTmp]), ncol= 5)
    colnames(tableOut) <- c('signature','Wilcox_pval','q25','q50','q75')
    tableOut[,1] <- paste(paste('c',contrast,sep=''),regulation[compTmp],sep='_')
  
    #Calculate percentiles for Background
    tmpBg <- sort(as.numeric(codonCalcOut))
    ecdfBg <- 1:length(tmpBg)/length(tmpBg)
    bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
    bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
    bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
  
  #
    for(i in 1:2){
      tmpVal <- as.numeric(codonCalcOut[names(codonCalcOut) %in% resComp[[i]]])
      lines(ecdf(tmpVal),col=as.character(AnotaColours[regulation[compTmp][i]]),main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=4)
    #
      tableOut[i,2] <- format(as.numeric(wilcox.test(tmpVal,as.numeric(codonCalcOut),alternative='two.sided')[3]),scientific = TRUE,digits=2)
      #Calculate percentiles for signatures and difference from background
      tmpSign <- sort(tmpVal)
      ecdfSign <- 1:length(tmpSign)/length(tmpSign)
      tableOut[i,3] <- format(tmpSign[which(ecdfSign >= 0.25)[1]] - bg_025, digits = 2)
      tableOut[i,4] <- format(tmpSign[which(ecdfSign >= 0.5)[1]] - bg_05, digits = 2)
      tableOut[i,5] <- format(tmpSign[which(ecdfSign >= 0.75)[1]] - bg_075, digits = 2)
    }
    plotrix::addtable2plot(0,1.01,tableOut,bty="n",display.rownames=FALSE,hlines=FALSE,vlines=TRUE,title="",cex = 1,bg=as.character(AnotaColours[regulation[compTmp]]),xpad=0.2,ypad=1.4)
    dev.off()
  }
  #
  return(codonCalcOut)
}