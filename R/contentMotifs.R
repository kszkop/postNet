##Run analysis length
contentMotifs <- function(ads=NULL,
                          regulation=NULL,
                          contrast=NULL,
                          comparisons=NULL,
                          customBg=NULL,
                          geneList=NULL,
                          geneListcolours=NULL,
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
                          pdfName=NULL,
                          plotOut=TRUE
){
  
  #Subset annot for only expressed genes
  annotBg <- gSel(annot=annot, ads=ads, customBg=customBg, geneList=geneList)
  #Select region of interest
  annotTmp <- regSel(annot=annotBg, region=region)
  #Per gene
  annotBgSel <- isoSel(annot=annotTmp, method=selection)
  
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
    motifOutTmp <- as.numeric(sapply(annotBgSel$seqTmp, calc_g4, min_score = min_score))
    names(motifOutTmp) <- as.character(annotBgSel$geneID)
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
  if(isTRUE(plotOut)){
    #
    resOut <- resSel(vIn=motifOut, ads=ads, regulation=regulation, contrast=contrast, customBg=customBg, geneList=geneList)
    #
    coloursOut <- coloursSel(ads=ads, regulation=regulation, resOut=resOut, geneList=geneList, geneListcolours=geneListcolours,customBg=customBg)
    #
    #Plot
    pdf(ifelse(is.null(pdfName),paste(region,motif,'content.pdf',sep='_'), paste(pdfName,region,motif,'content.pdf',sep='_')),width= 8,height=8, useDingbats = F)
    par(mar=c(5,5,8,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
    #
    xlim_min <- ifelse(isTRUE(resid),floor(quantile(as.numeric(unlist(resOut)),0.01)), 0)
    xlim_max <- roundUpNice(abs(as.numeric(quantile(as.numeric(unlist(resOut)),0.99))))
        
    plot(ecdf(resOut[[1]]),col=coloursOut[1],main='',xlab='',ylab='',verticals=TRUE, do.p=FALSE,lwd=3,bty="n",yaxt="none",font=2, xlim=c(xlim_min,xlim_max),xaxt="none")
    #
    for(i in 2:length(resOut)){
      lines(ecdf(resOut[[i]]),col=coloursOut[i],main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=4)
    }
    
    mtext(side=1, line=4, paste('Number of motifs \n',motif,sep=''), col='black', font=2,cex=1.2)
    mtext(side=2, line=3, 'Fn(x)', col="black", font=2, cex=1.2)
    
    axis(side=1,seq(xlim_min ,xlim_max,1),font=2,lwd=2)
    axis(side=2,seq(0,1,0.2), font=2,las=2,lwd=2)
    #
    if(!is.null(comparisons)){
      #Prepare comparisons per regulation
      for(j in 1:length(comparisons)){
        if(!is.null(ads) | !is.null(customBg)){
          compTmp <- comparisons[[j]]+1
        } else {
          compTmp <- comparisons[[j]]
        }
        #stats
        pvalTmp <- format(as.numeric(wilcox.test(resOut[[compTmp[1]]], resOut[[compTmp[2]]],alternative='two.sided')[3]),scientific = TRUE,digits=2)
        #
        tableOut <- matrix(NA, nrow=length(comparisons), ncol= 5)
        colnames(tableOut) <- c('signature','Wilcox_pval','q25','q50','q75')
        tableOut[j,1] <- paste(names(resOut)[compTmp[2]],'vs', names(resOut)[compTmp[1]],sep=' ')
        tableOut[j,2] <- pvalTmp
        
        #Calculate percentiles
        tmpBg <- sort(resOut[[compTmp[1]]])
        ecdfBg <- 1:length(tmpBg)/length(tmpBg)
        bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
        bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
        bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
        
        #Calculate percentiles for second and difference from background
        tmpSign <- sort(resOut[[compTmp[2]]])
        ecdfSign <- 1:length(tmpSign)/length(tmpSign)
        tableOut[j,3] <- format(tmpSign[which(ecdfSign >= 0.25)[1]] - bg_025, digits = 2)
        tableOut[j,4] <- format(tmpSign[which(ecdfSign >= 0.5)[1]] - bg_05, digits = 2)
        tableOut[j,5] <- format(tmpSign[which(ecdfSign >= 0.75)[1]] - bg_075, digits = 2)
        #
        if(length(which(grepl('background',c(names(resOut)[compTmp[2]], names(resOut)[compTmp[1]]))))>0){
          colT <- gsub("\\_.*","", names(resOut)[compTmp][which(names(resOut)[compTmp] != 'background')])
          colT <- coloursOut[colT]
        } else {
          colT <- 'white'
        }
        plotrix::addtable2plot(xlim_min,1.01,tableOut,bty="n",display.rownames=FALSE,hlines=FALSE,vlines=TRUE,title="",cex = 0.9,bg=colT,xpad=0.1,ypad=1.4,xjust=0,yjust=1)
      }
    }
    dev.off()
    #
  }
  return(motifOut)
}
  
  
  