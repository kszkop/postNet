signaturesHeatmap <- function(ptn,
                              signatureList,
                              unit = 'FDR',
                              addEff = NULL,
                              pdfName = NULL){
  #
  check_ptn(ptn)
  check_geneList(signatureList)
  if(!check_shiftUnit(unit)){
    stop("'unit' can be only 'FDR' or p25, p50, p75 or any percentile in format p with number). 
         If FDR: -log10 FDR wilcoxon test corrected for multitesting with direction of ecdf shift, 
         if percentile: difference from background for one of the percentiles from ecdf")
  }
  if(!is.null(addEff)){
    if(!is_named_list_of_named_numeric_vectors(addEff)){
      stop('addEff must be a named list of named numeric vectors')
    }
  }
  compOut <- ifelse(is.null(addEff), 1, length(addEff)+1)
  effNames <- ifelse(is.null(addEff), 'ptn_effect', c('ptn_effect', names(addEff)))
  #
  tableFinal <- matrix(NA, nrow=length(signatureList), ncol=compOut)
  row.names(tableFinal) <- names(signatureList)
  colnames(tableFinal) <- effNames

  for(i in 1:compOut){
    if(i == 1){
      effIn <- ptn_effect(ptn)
    } else {
      effIn <- addEff[[i]]
    }
    regData <- data.frame(geneSymb = names(effIn))
    regData$effIn <- as.numeric(effIn)
    
    ##colculate metric for each signature
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
        pval <- as.numeric(wilcox.test(as.numeric(regData[regData$signature==names(signatureList)[sign],]$effIn),as.numeric(regData[regData$signature=='bkg',]$effIn),alternative='two.sided')[3])
        if(pval == 0){
          pval <- 1e-300
        }
        fdrOut[sign] <- pval
      }
    }
    if(unit=='FDR'){
      adjFDR <- -log10(p.adjust(fdrOut))
      fdrDirec <- ifelse(percOut>=0, adjFDR*1,adjFDR*-1)
      tableFinal[,i] <- fdrDirec
    } else {
      tableFinal[,i] <- percOut
    }
  }
  if(unit=='FDR'){
    if(max(abs(c(min(as.vector(tableFinal)),max(as.vector(tableFinal)))))>10){
      breaks <- seq(-10,10,length.out=20)
    } else {
      breaks <- seq(min(as.vector(tableFinal)),max(as.vector(tableFinal)), length.out=25)
    }
    keyL <- paste('-log10',unit,sep=' ')
  } else {
    breaks <- seq(-max(abs(min(as.vector(tableFinal))),abs(max(as.vector(tableFinal)))),max(abs(min(as.vector(tableFinal))),abs(max(as.vector(tableFinal)))), length.out=25)
    keyL <- unit
  }
  
  breaks <- sort(unique(c(breaks,0)))
  len <- length(breaks) -1 
  
  if(length(compOut)==1){
    tableFinal <- cbind(tableFinal,rep(0,nrow(tableFinal)))
  }
  
  colTmp <- rev(colorRampPalette(RColorBrewer::brewer.pal(name="RdBu",n=11))(len))
  colTmp[(len + 1) / 2] <- "#FFFFFF"
  
  pdf(ifelse(is.null(pdfName),'heatmap.pdf', paste(pdfName,'heatmap.pdf',sep='_')),width= 8,height=8, useDingbats = F)
  par(mar=c(10,5,5,10),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.9, cex.main=0.7,cex.lab=0.9)

  gplots::heatmap.2(tableFinal,trace = "none",
                    breaks = breaks,
                    Rowv = TRUE,
                    col = colTmp,
                    key.xlab = keyL,key.title="",dendrogram   = "none", Colv=FALSE, ,density.info = "none",
                    tracecol = NULL,margins = c(10,19),lhei=c(1,6),lwid=c(0.5,1),
                    cexRow = 0.9,cexCol = 0.9,offsetRow = -5)
  dev.off()
}
