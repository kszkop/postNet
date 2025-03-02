plotSignatures <- function(ptn,
                           signatureList,
                           signature_colours=NULL, 
                           xlim = c(-2,2), 
                           tableCex = 1 ,
                           pdfName = NULL){
  #
  check_ptn(ptn)
  check_geneList(signatureList)

  signNames <- names(signatureList)
  if(!is.null(signature_colours)){
    if (!is.character(signature_colours) || !length(signature_colours)== length(signatureList)) {
      stop("'signature_colours' should be a character vector of the same length as number of signatures in signatureList. These colours will be used for plotting.")
    }
  } else {
    signature_colours <- paste0("#", sprintf("%06X", sample(0:16777215, length(signatureList))))
  }

  if(!check_number(tableCex)){
    stop("please provide number for tableCex to scale size of table with statistics")
  }

  ##apvEff of effect
  effIn <- ptn_effect(ptn)
  regData <- data.frame(geneSymb = names(effIn))
  regData$effIn <- as.numeric(effIn)
  #regData <- data.frame(geneSymb = rownames(ads@dataP))
  
  #regData$totalApvEff <- ads@totalmRNA@apvStatsRvm[[contrast]][,"apvEff"]
  #regData$polyApvEff <- ads@translatedmRNA@apvStatsRvm[[contrast]][,"apvEff"]
  #regData$buffApvEff <- ads@buffering@apvStatsRvm[[contrast]][,"apvEff"]
  #regData$translationApvEff <- ads@translation@apvStatsRvm[[contrast]][,"apvEff"]
  
  #Select for scatter
  if(is.null(scatterXY)){
    scatterXY <- roundNice(max(abs(c(range(regData$totalApvEff),range(regData$polyApvEff)))),direction='up')
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
        xmin <- ifelse(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.01))<0, -roundNice(abs(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.01))),direction='up'),roundNice(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.01))),direction='up')
        xmax <- ifelse(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.99))<0, -roundNice(abs(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.01))),direction='up'),roundNice(as.numeric(quantile(as.numeric(unlist(tmpBgOut)),0.99))),direction='up')
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
        xmin <- ifelse(as.numeric(quantile(as.numeric(tmpBg),0.01))<0, -roundNice(abs(as.numeric(quantile(as.numeric(tmpBg),0.01))),direction='up'),roundNice(as.numeric(quantile(as.numeric(tmpBg),0.01))),direction='up')
        xmax <- ifelse(as.numeric(quantile(as.numeric(tmpBg),0.99))<0, -roundNice(abs(as.numeric(quantile(as.numeric(tmpBg),0.01))),direction='up'),roundNice(as.numeric(quantile(as.numeric(tmpBg),0.99))),direction='up')
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