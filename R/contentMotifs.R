##Run analysis length
contentMotifs <- function(ads,
                            regulation,
                            contrast,
                            comparisons=NULL,
                            motif, #if want G4 write G4
                            seqType='dna', #protein/dna/rna so is matching motif
                            len=1,
                            min_score=47, 
                            region, #UTR5, CDS, UTR3
                            subregion=NULL, #number of nucleotides from start if positive or end if negative.
                            subregionSel, #select or exclude , required if subregion is not null.
                            annot, 
                            selection, #shortest, longest, random (default)
                            resid=FALSE,
                            pdfName=NULL
){
  #Subset annot for only expressed genes
  bg <- row.names(ads@dataP)
  annotBg <- annot[annot$geneID %in% bg,]
  
  #Select region of interest
  if(region=='UTR5'){
    seqTmp <- annotBg$UTR5_seq
    lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
    annotBg <- cbind(annotBg[,c(1:2)], seqTmp,lenTmp)
  }
  if(region=='UTR3'){
    seqTmp <- annotBg$UTR3_seq
    lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
    annotBg <- cbind(annotBg[,c(1:2)], seqTmp,lenTmp)
  }
  if(region=='CDS'){
    seqTmp <- annotBg$CDS_seq
    lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
    annotBg <- cbind(annotBg[,c(1:2)], seqTmp,lenTmp)
  }
  
  #Select per gene level
  if(selection=='shortest'){
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice(which.min(lenTmp)))
  } else if(selection=='longest'){
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice(which.max(lenTmp)))
  } else {
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice_sample(n = 1))
  }
  #
  if(seqType=='protein'){
    proseqtmp <- as.character(sapply(annotBgSel$seqTmp, function(x) seqinr::c2s(seqinr::translate(seqinr::s2c(x)))))
    #
    annotBgSel$seqTmp <- proseqtmp
    #
    annotBgSel$lenTmp <- annotBgSel$lenTmp/3
  }
  #
  if(!is.null(subregion)){
    #
    subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos=subregion,subregionSel=subregionSel)))
    #
    annotBgSel$seqTmp <- subSeq
  }
  annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp),]
  
  #Find given motif
  if(motif=="G4"){
    motifOut <- as.numeric(sapply(annotBgSel$seqTmp, calc_g4, min_score = min_score))
    names(motifOut) <- as.character(annotBgSel$geneID)
  } else {
    motif <- toupper(motif)
    if(seqType=='dna' | seqType=='rna'){
      motifTmp <- convertIUPAC(motif)
    } else {
      motifTmp <- replaceProtAmbig(motif)
    }
    #calcuate distance
    len <- ifelse(len==1, length(seqinr::s2c(motifTmp)), (length(seqinr::s2c(motifTmp)))+(len-1))
    #
    motifOutTmp <- as.numeric(sapply(annotBgSel$seqTmp, calc_motif, motif=motifTmp, len=len))
    names(motifOutTmp) <- as.character(annotBgSel$geneID)
  }

  #
  if(isTRUE(resid)){
    motifOut <- lm(as.numeric(motifOutTmp) ~ log2(as.numeric(annotBgSel$lenTmp)))$residuals
    names(motifOut) <- names(motifOutTmp)
  } else {
    motifOut <- motifOutTmp
  }
  
  
  #
  #Extract all results
  results <- anota2seqGetDirectedRegulations(ads)
  
  #Prepare plotting
  AnotaColours <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],RColorBrewer::brewer.pal(8,"Reds")[c(2,6)],RColorBrewer::brewer.pal(8,"Greens")[c(4,8)], RColorBrewer::brewer.pal(8,"Greens")[c(2,6)],RColorBrewer::brewer.pal(8,"Blues")[c(4,8)])
  names(AnotaColours) <- c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown","bufferingmRNAUp","bufferingmRNADown")
  #subset for desired regulations depending on number of contrasts
  #if(length(unique(contrast))==1){
  #  res <- results[[1]][names(results[[1]]) %in% regulation]
  #} else {
  res <- vector("list", length = length(regulation))
  for(i in unique(contrast)){
    resTmp <- results[[i]][names(results[[i]]) %in% regulation[which(contrast==i)]]
    res[which(contrast==i)] <- resTmp
  }
  #
  #Prepare comparisons per regulation
  for(j in 1:length(comparisons)){
    compTmp <- comparisons[[j]]
    #
    resComp <- res[compTmp]
    #Plot
    pdf(ifelse(is.null(pdfName),paste(region,motif,paste(regulation[compTmp],collapse='_'),'content.pdf',sep='_'), paste(pdfName,region,motif,paste(regulation[compTmp],collapse='_'),'content.pdf',sep='_')),width= 8,height=8, useDingbats = F)
    par(mar=c(5,5,8,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
    plot(ecdf(as.numeric(motifOut)),col='grey45',main='',xlab='',ylab='',verticals=TRUE, do.p=FALSE,lwd=3,bty="n",yaxt="none",font=2, xlim=c(0,as.numeric(quantile(as.numeric(motifOut),0.99))),xaxt="none")
  
    mtext(side=1, line=4, paste('Number of motifs \n',motif,sep=''), col='black', font=2,cex=1.2)
    mtext(side=2, line=3, 'Fn(x)', col="black", font=2, cex=1.2)
  
    axis(side=1,seq(0,as.numeric(quantile(as.numeric(motifOut),0.99)),1), font=2,lwd=2)
    axis(side=2,seq(0,1,0.2), font=2,las=2,lwd=2)
  
    #
    tableOut <- matrix(NA, nrow= length(regulation[compTmp]), ncol= 5)
    colnames(tableOut) <- c('signature','Wilcox_pval','q25','q50','q75')
    tableOut[,1] <- paste(paste('c',contrast,sep=''),regulation[compTmp],sep='_')
  
    #Calculate percentiles for Background
    tmpBg <- sort(as.numeric(motifOut))
    ecdfBg <- 1:length(tmpBg)/length(tmpBg)
    bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
    bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
    bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
  
    #
    for(i in 1:2){
      tmpVal <- as.numeric(motifOut[names(motifOut) %in% resComp[[i]]])
      lines(ecdf(tmpVal),col=as.character(AnotaColours[regulation[compTmp][i]]),main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=4)
    #
      tableOut[i,2] <- format(as.numeric(wilcox.test(tmpVal,as.numeric(motifOut),alternative='two.sided')[3]),scientific = TRUE,digits=2)
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
  return(motifOut) 
}
  
  
  
  