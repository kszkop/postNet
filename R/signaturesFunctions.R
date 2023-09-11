signatureFunction <- function(signatureList,
                              generalName,
                              dataName,
                              effects_names=c('total mRNA log2FC', 'polysome associated mRNA log2FC','buffering log2FC','translation log2FC'),
                              colours, 
                              ads,
                              contrast,
                              xlim=NULL, 
                              scatterXY=NULL,
                              tableCex,
                              pdfName){
  #
  signNames <- names(signatureList)
  ##apvEff of effect
  regData <- data.frame(geneSymb = rownames(ads@dataP))
  
  regData$totalApvEff <- ads@totalmRNA@apvStatsRvm[[contrast]][,"apvEff"]
  regData$polyApvEff <- ads@translatedmRNA@apvStatsRvm[[contrast]][,"apvEff"]
  regData$buffApvEff <- ads@buffering@apvStatsRvm[[contrast]][,"apvEff"]
  regData$translationApvEff <- ads@translation@apvStatsRvm[[contrast]][,"apvEff"]
  
  #Select for scatter
  if(is.null(scatterXY)){
    scatterXY <- roundUpNice(max(abs(c(range(regData$totalApvEff),range(regData$polyApvEff)))))
  }
  #
  pdf(ifelse(is.null(pdfName),paste('data_',dataName,'_signature_',generalName,'.pdf',sep=''), paste(pdfName,paste('data_',dataName,'_signature_',generalName,'.pdf',sep=''),sep='_')),width= 18,height=4, useDingbats = F)
  par(mfrow=c(1,5),mar=c(5,5,6,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.25,cex.main=1.9,cex.lab=2)
  plot(regData$totalApvEff,regData$polyApvEff,pch=16,cex=1.9,col='grey75',ylab=effects_names[2],xlab=effects_names[1], main=paste(paste('Data: ',dataName,sep=''),'\n',paste('Signature: ',generalName,''),sep=''), xlim=c(-scatterXY,scatterXY), ylim=c(-scatterXY,scatterXY))
  abline(v=0)
  abline(h=0)
  
  #Check whether gene list overlap
  if(any(duplicated(unlist(signatureList)))){
    cat("There are some genes that overlap between signatures. Separate background for each will be used")
    ##collect signatures
    for(i in 1:length(signatureList)){
      regData[,5+i] <- 'bkg'
      regData[,5+i][regData$geneSymb %in% signatureList[[i]]] <- names(signatureList)[i]
      #plot 
      points(regData$totalApvEff[regData[,5+i]==names(signatureList)[i]],regData$polyApvEff[regData[,5+i]==names(signatureList)[i]],col=colours[i],pch=16,cex=1.7)
    }
  } else {
    regData$signature <- 'bkg'
    for(i in 1:length(signatureList)){
      regData$signature[regData$geneSymb %in% signatureList[[i]]] <- names(signatureList)[i]
      #plot 
      points(regData$totalApvEff[regData$signature==names(signatureList)[i]],regData$polyApvEff[regData$signature==names(signatureList)[i]],col=colours[i],pch=16,cex=1.7)
    }
  }
  #
  #if(!is.null(signatureListAdd)){
  #  regData$signatureAdd  <- 'bkg'
  #}
  #if(!is.null(signatureListAdd)){
  #  for(i in 1:length(signatureListAdd)){
  #    regData$signatureAdd[regData$geneSymb %in% signatureListAdd[[i]]] <- signatureListAddNames[i]
  #  }
  #}
  
  for(eff in 5:2){
    ##Calculate statistics
    tableOut <- matrix(NA, nrow= length(signNames), ncol= 5)
    colnames(tableOut) <- c('signature','Wilcox_pval','q25','q50','q75')
    tableOut[,1] <- signNames
    
    if(any(duplicated(unlist(signatureList)))){
      tmpBgOut <- list()
      for(i in 1:length(signatureList)){
        #Calculate percentiles for Background
        tmpBg <- sort(as.numeric(regData[regData[,(5+i)]=='bkg',][,eff]))
        ecdfBg <- 1:length(tmpBg)/length(tmpBg)
        bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
        bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
        bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
        #
        tmpBgOut[[i]] <- tmpBg
        #
        tableOut[i,2] <- format(as.numeric(wilcox.test(as.numeric(regData[regData[,(5+i)]==signNames[i],][,eff]),as.numeric(regData[regData[,(5+i)]=='bkg',][,eff]),alternative='two.sided')[3]),scientific = TRUE,digits=2)
        
        #Calculate percentiles for signatures and difference from background
        tmpSign <- sort(as.numeric(regData[regData[,(5+i)]==signNames[i],][,eff]))
        ecdfSign <- 1:length(tmpSign)/length(tmpSign)
        
        tableOut[i,3] <- format(tmpSign[which(ecdfSign >= 0.25)[1]] - bg_025, digits = 2)
        tableOut[i,4] <- format(tmpSign[which(ecdfSign >= 0.5)[1]] - bg_05, digits = 2)
        tableOut[i,5] <- format(tmpSign[which(ecdfSign >= 0.75)[1]] - bg_075, digits = 2)
        #
      }
      #Plot ecdfs
      if(!is.null(xlim)){
        xmin <- xlim[1]
        xmax <- xlim[2]
      } else {
        xmin <- ifelse(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.01))<0, -roundUpNice(abs(as.numeric(quantile(as.numeric(tmpBg),0.01)))),roundUpNice(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.01))))
        xmax <- ifelse(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.99))<0, -roundUpNice(abs(as.numeric(quantile(as.numeric(tmpBg),0.01)))),roundUpNice(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.99))))
      }
      plot(ecdf(as.numeric(regData[regData[,(5+i)]=='bkg',][,eff])),col='grey55',main='',xlab=effects_names[eff-1],verticals=TRUE, do.p=FALSE,lwd=3,xlim=c(xmin,xmax))
      legend(xmin,0.95,fill='grey55',border='grey55','Background',bty='n',cex=1.3)
      
      #
      for(i in 1:length(signatureList)){
        lines(ecdf(as.numeric(regData[regData[,(5+i)]==signNames[i],][,eff])),col=colours[i],main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=3)
      }
      plotrix::addtable2plot(xmin-abs((xmin*0.1)),1.05,tableOut,bty="n",display.rownames=FALSE,hlines=FALSE,vlines=TRUE,title="",cex = tableCex,bg=colours,xpad=0.2,ypad=1.4)
    } else {
      #Calculate percentiles for Background
      tmpBg <- sort(as.numeric(regData[regData$signature=='bkg',][,eff]))
      ecdfBg <- 1:length(tmpBg)/length(tmpBg)
      bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
      bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
      bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
      
      #Calculate Wilcoxon pvalue
      for(i in 1:length(signatureList)){
        tableOut[i,2] <- format(as.numeric(wilcox.test(as.numeric(regData[regData$signature==signNames[i],][,eff]),as.numeric(regData[regData$signature=='bkg',][,eff]),alternative='two.sided')[3]),scientific = TRUE,digits=2)
      
        #Calculate percentiles for signatures and difference from background
        tmpSign <- sort(as.numeric(regData[regData$signature==signNames[i],][,eff]))
        ecdfSign <- 1:length(tmpSign)/length(tmpSign)
      
        tableOut[i,3] <- format(tmpSign[which(ecdfSign >= 0.25)[1]] - bg_025, digits = 2)
        tableOut[i,4] <- format(tmpSign[which(ecdfSign >= 0.5)[1]] - bg_05, digits = 2)
        tableOut[i,5] <- format(tmpSign[which(ecdfSign >= 0.75)[1]] - bg_075, digits = 2)
      }
      #Plot ecdfs
      if(!is.null(xlim)){
        xmin <- xlim[1]
        xmax <- xlim[2]
      } else {
        xmin <- ifelse(as.numeric(quantile(as.numeric(tmpBg),0.01))<0, -roundUpNice(abs(as.numeric(quantile(as.numeric(tmpBg),0.01)))),roundUpNice(as.numeric(quantile(as.numeric(tmpBg),0.01))))
        xmax <- ifelse(as.numeric(quantile(as.numeric(tmpBg),0.99))<0, -roundUpNice(abs(as.numeric(quantile(as.numeric(tmpBg),0.01)))),roundUpNice(as.numeric(quantile(as.numeric(tmpBg),0.99))))
      }
      plot(ecdf(as.numeric(regData[regData$signature=='bkg',][,eff])),col='grey55',main='',xlab=effects_names[eff-1],verticals=TRUE, do.p=FALSE,lwd=3,xlim=c(xmin,xmax))
      legend(xmin,0.95,fill='grey55',border='grey55','Background',bty='n',cex=1.3)
      #
      for(i in 1:length(signatureList)){
        lines(ecdf(as.numeric(regData[regData$signature==signNames[i],][,eff])),col=colours[i],main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=3)
      }
      plotrix::addtable2plot(xmin-abs((xmin*0.1)),1.05,tableOut,bty="n",display.rownames=FALSE,hlines=FALSE,vlines=TRUE,title="",cex = tableCex,bg=colours,xpad=0.2,ypad=1.4)
    }
    #if(!is.null(signatureListAddNames)){
    #  legend(xmin,0.9,fill=coloursAdd,border=coloursAdd,'Regulated',bty='n',cex=1.3)
    #}
    #if(!is.null(signatureListAdd)){
    #  for(i in 1:length(signatureListAdd)){
    #    lines(ecdf(as.numeric(regData[regData$signatureAdd==signatureListAddNames[i],][,eff])),col=coloursAdd[i],main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=1,lty='dashed')
    #  }
    #}
  }
  dev.off()
}

###
signaturesHeatmap <- function(signatureList, #allSignatures
                                      ads, #output from anota2seq
                                      contrast, #it can be vector of contrasts
                                      contrastNames, #it can be also a vector of equal length to contrasts
                                      unit, #(options: FDR, p25, p50, p75 or any percentile in format p with number), what will be used: -log10 FDR wilcoxon test corrected for multitesting with direction of ecdf shift, or difference from background for one of the percentiles from ecdf.
                                      RegMode,#(options: total, poly, translation), whether respectively total mRNA log2FC, poly mRNA log2FC or translation log2 FC
                                      pdfName
){
  #For each contrast:
  tableFinal <- matrix(NA, nrow=length(signatureList), ncol=length(contrasts))
  row.names(tableFinal) <- names(signatureList)
  colnames(tableFinal) <- contrastNames
  
  for(cont in 1:length(contrast)){
    ##apvEff of effect
    regData <- data.frame(geneSymb = rownames(ads@dataP))
    if(RegMode=='translation'){
      regData$ApvEff<- ads@translation@apvStatsRvm[[cont]][,"apvEff"]
    }
    if(RegMode=='buffering'){
      regData$ApvEff<- ads@buffering@apvStatsRvm[[cont]][,"apvEff"]
    }
    if(RegMode=='poly'){
      regData$ApvEff <- ads@translatedmRNA@apvStatsRvm[[cont]][,"apvEff"]
    }
    if(RegMode=='total'){
      regData$ApvEff <- ads@totalmRNA@apvStatsRvm[[cont]][,"apvEff"]
    }
    
    ##colculate metric for each signature
    percOut <- as.numeric()
    if(unit=='FDR'){
      fdrOut <- as.numeric()
    }
    for(sign in 1:length(signatureList)){
      regData[,3] <- 'bkg'
      regData[,3][regData$geneSymb %in% signatureList[[sign]]] <- names(signatureList)[sign]
      colnames(regData)[3] <- 'signature'
      
      #Calculate percentiles for Background
      tmpBg <- sort(as.numeric(regData[regData$signature=='bkg',]$ApvEff))
      ecdfBg <- 1:length(tmpBg)/length(tmpBg)
      if(unit=='FDR'){
        percentileBG <- tmpBg[which(ecdfBg >= 0.5)[1]]
      } else {
        percentileBG <- tmpBg[which(ecdfBg >= as.numeric(gsub('p','0.',unit)))[1]]
      }
      tmpSign <- sort(as.numeric(regData[regData$signature==names(signatureList)[sign],]$ApvEff))
      ecdfSign <- 1:length(tmpSign)/length(tmpSign)
      if(unit=='FDR'){
        percentileSign <- tmpSign[which(ecdfSign >= 0.5)[1]]
      } else {
        percentileSign <- tmpSign[which(ecdfSign >= as.numeric(gsub('p','0.',unit)))[1]]
      }
      percentileDiff <- percentileSign - percentileBG
      percOut[sign] <- percentileDiff
      #If metric chosen is FDR
      if(unit=='FDR'){
        #wilcox test
        pval <- as.numeric(wilcox.test(as.numeric(regData[regData$signature==names(signatureList)[sign],]$ApvEff),as.numeric(regData[regData$signature=='bkg',]$ApvEff),alternative='two.sided')[3])
        fdrOut[sign] <- pval
      }
    }
    if(unit=='FDR'){
      adjFDR <- -log10(p.adjust(fdrOut))
      fdrDirec <- ifelse(percOut>=0, adjFDR*1,adjFDR*-1)
      tableFinal[,cont] <- fdrDirec
    } else {
      tableFinal[,cont] <- percOut
    }
  }
  if(unit=='FDR'){
    if(max(abs(c(min(as.vector(tableFinal)),max(as.vector(tableFinal)))))>10){
      breaks <- seq(-10,10,length.out=20)
      len <- length(breaks)-1
    } else {
      breaks <- seq(min(as.vector(tableFinal)),max(as.vector(tableFinal)), length.out=25)
      len <- length(breaks)-1
    }
    keyL <- paste('-log10',unit,sep=' ')
  } else {
    breaks <- seq(-max(abs(min(as.vector(tableFinal))),abs(max(as.vector(tableFinal)))),max(abs(min(as.vector(tableFinal))),abs(max(as.vector(tableFinal)))), length.out=25)
    len <- length(breaks)-1
    keyL <- unit
  }
  if(length(contrast)==1){
    tableFinal <- cbind(tableFinal,rep(0,nrow(tableFinal)))
  }
  
  pdf(ifelse(is.null(pdfName),'heatmap.pdf', paste(pdfName,'heatmap.pdf',sep='_')),width= 8,height=8, useDingbats = F)
  par(mar=c(10,5,5,10),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.9, cex.main=0.7,cex.lab=0.9)
  gplots::heatmap.2(tableFinal,trace = "none",
            breaks = breaks,
            Rowv = TRUE,
            col = rev(colorRampPalette(RColorBrewer::brewer.pal(name="RdBu",n=11))(len)),
            key.xlab = keyL,key.title="",dendrogram   = "row",density.info = "none",
            tracecol = NULL,margins = c(10,19),lhei=c(1,6),lwid=c(0.5,1),
            cexRow = 0.9,cexCol = 0.9)
  dev.off()
}
