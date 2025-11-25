signaturesHeatmap <- function(ptn,
                              signatureList,
                              unit = 'FDR',
                              pdfName = NULL){
  #
  check_ptn(ptn)
  check_geneList(signatureList)
  if(!check_shiftUnit(unit)){
    stop("'unit' can be only 'FDR' or p25, p50, p75 or any percentile in format p with number). 
         If FDR: -log10 FDR wilcoxon test corrected for multitesting with direction of ecdf shift, 
         if percentile: difference from background for one of the percentiles from ecdf")
  }
  
  compOut <- 1
  effNames <- 'ptn_effect'
  #
  tableFinal <- matrix(NA, nrow=length(signatureList), ncol=compOut)
  row.names(tableFinal) <- names(signatureList)
  colnames(tableFinal) <- effNames

  
  effIn <- ptn_effect(ptn)
  regData <- data.frame(geneSymb = names(effIn))
  regData$effIn <- as.numeric(effIn)
    
  ##calculate metric for each signature
  percOut <- as.numeric()
  if(unit=='FDR'){
    fdrOut <- as.numeric()
  }
  for(sign in 1:length(signatureList)){
    regData[,3] <- 'bkg'
    regData[,3][regData$geneSymb %in% signatureList[[sign]]] <- names(signatureList)[sign]
    colnames(regData)[3] <- 'signature'
      
    #
    tmpBg <- sort(as.numeric(regData[regData$signature=='bkg',]$effIn))
    ecdfBg <- 1:length(tmpBg)/length(tmpBg)
    if(unit=='FDR'){
      percentileBG <- tmpBg[which(ecdfBg >= 0.5)[1]]
    } else {
      percentileBG <- tmpBg[which(ecdfBg >= as.numeric(gsub('p','0.',unit)))[1]]
    }
    tmpSign <- sort(as.numeric(regData[regData$signature==names(signatureList)[sign],]$effIn))
    ecdfSign <- 1:length(tmpSign)/length(tmpSign)
    if(unit=='FDR'){
      percentileSign <- tmpSign[which(ecdfSign >= 0.5)[1]]
    } else {
      percentileSign <- tmpSign[which(ecdfSign >= as.numeric(gsub('p','0.',unit)))[1]]
    }
    percentileDiff <- percentileSign - percentileBG
    percOut[sign] <- percentileDiff

    if(unit=='FDR'){
      #wilcox test
      pval <- as.numeric(wilcox.test(as.numeric(regData[regData$signature==names(signatureList)[sign],]$effIn),as.numeric(regData[regData$signature=='bkg',]$effIn),exact=FALSE, alternative='two.sided')[3])
      if(pval == 0){
        pval <- 1e-300
      }
      fdrOut[sign] <- pval
    }
  }
  if(unit=='FDR'){
    adjFDR <- -log10(p.adjust(fdrOut))
    fdrDirec <- ifelse(percOut>=0, adjFDR*1,adjFDR*-1)
    tableFinal[,1] <- fdrDirec
  } else {
    tableFinal[,1] <- percOut
  }

  breaks <- seq(-max(abs(min(as.vector(tableFinal))),abs(max(as.vector(tableFinal)))),max(abs(min(as.vector(tableFinal))),abs(max(as.vector(tableFinal)))), length.out=51)
  keyL <- ifelse(unit=='FDR', paste('-log10',unit,sep=' '), unit)

  colTmp <- grDevices::colorRampPalette(c("blue", "white", "red"))(50)
    
  if(length(compOut)==1){
    tableFinal <- cbind(tableFinal,rep(0,nrow(tableFinal)))
  }

  pdf(ifelse(is.null(pdfName),'heatmap.pdf', paste(pdfName,'heatmap.pdf',sep='_')),width= 8,height=8, useDingbats = F)
  par(mar=c(10,5,5,10),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.9, cex.main=0.7,cex.lab=0.9)

  gplots::heatmap.2(tableFinal,trace = "none",
                    breaks = breaks,
                    Rowv = TRUE,
                    col = colTmp,
                    key.xlab = keyL,key.title="",dendrogram   = "row", Colv=FALSE, ,density.info = "none",
                    tracecol = NULL,margins = c(10,19),lhei=c(1,6),lwid=c(0.5,1),
                    cexRow = 0.9,cexCol = 0.9,offsetRow = -5)
  dev.off()
}
