##Run analysis length
lengthAnalysis <- function(ads,
                           regulation,
                           contrast,
                           comparisons,
                           region, #UTR5, CDS, UTR3
                           annot,
                           selection, #shortest, longest, random (default)
                           plotType='boxplot',# option 'violin' or 'ecdf'
                           geneList=NULL,
                           geneListnames=NULL,
                           geneListcolours=NULL,
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
  } else if (selection=='longest'){
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice(which.max(lenTmp)))
  } else {
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice_sample(n = 1))
  }
  #
  lenForAnalysis <- annotBgSel$lenTmp
  names(lenForAnalysis) <- annotBgSel$geneID
  
  #Prepare plotting
  AnotaColours <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],RColorBrewer::brewer.pal(8,"Greens")[c(4,8)], RColorBrewer::brewer.pal(8,"Blues")[c(4,8)])
  names(AnotaColours) <- c("translationUp","translationDown", "mRNAAbundanceUp","mRNAAbundanceDown","bufferingmRNAUp","bufferingmRNADown")
  #
  #Extract all results
  results <- anota2seqGetDirectedRegulations(ads)
  
  #subset for desired regulations depending on number of contrasts
  #if(length(unique(contrast))==1){
  #  res <- results[[1]][names(results[[1]]) %in% regulation]
  #} else {
  res <- vector("list", length = length(regulation))
  for(i in unique(contrast)){
    resTmp <- results[[i]][names(results[[i]]) %in% regulation[which(contrast==i)]]
    res[which(contrast==i)] <- resTmp
  }
  #}
  if(!is.null(geneList)){
    res <- append(res, geneList)
  }
  #
  resCol <- list()
  #Load background first
  resCol[[1]] <- log2(as.numeric(lenForAnalysis))
  
  #Plot
  pdf(ifelse(is.null(pdfName),paste(region,plotType,'lengthAnalysis.pdf',sep='_'), paste(pdfName,region,plotType,'lengthAnalysis.pdf',sep='_')),width= 8,height=8, useDingbats = F)
  if(plotType=='boxplot'| plotType=='violin'){
    par(mar=c(8,12,5,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
    plot(1, 1, xlim=c(0.5,length(regulation)+ifelse(!is.null(geneList),length(geneList),0)+1.5), ylim=c(0,range(log2(as.numeric(lenForAnalysis)))[2]+(1.25*length(comparisons))), xaxt='n', xlab='', ylab='',type='n',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE)
  
    #axis(side=2,seq(0,round(range(log2(as.numeric(lenForAnalysis)))[2]+(1.25*length(comparisons))),2), font=2,las=2,lwd=2)
    axis(side=2, font=2,las=2,lwd=2,at=sapply(c(1,25,100,200,400,1000,4000,25000),log2),labels = c(0,25,100,200,400,1000,4000,25000))
  
    mtext(side=2, line=6, paste(region, 'Log2 length', sep=' '), col="black", font=2, cex=1.7,at=median(log2(as.numeric(lenForAnalysis))))
    if(!is.null(geneList)){
      text(1:(length(regulation)+length(geneList)+1), par("usr")[3] - 0.45, labels=c('background',paste(paste('c',contrast,sep=''),regulation,sep='_'), geneListnames), xpd=NA,cex=0.9,srt=45,adj=1)
    } else {
      text(1:(length(regulation)+1), par("usr")[3] - 0.45, labels=c('background',paste(paste('c',contrast,sep=''),regulation,sep='_')), xpd=NA,cex=0.9,srt=45,adj=1)
    }
    
    text(1,0,round(mean(as.numeric(lenForAnalysis)),0),font=2)
    abline(lty=5, h=median(log2(as.numeric(lenForAnalysis))))
    
  } else if (plotType=='ecdf'){
    #
    xlim_min <- as.numeric(quantile(log2(as.numeric(lenForAnalysis)),0.01))
    xlim_max <- as.numeric(quantile(log2(as.numeric(lenForAnalysis)),0.99))
    #
    par(mar=c(5,5,8,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
    plot(ecdf(log2(as.numeric(lenForAnalysis))),col='grey45',main='',xlab='',ylab='',verticals=TRUE, do.p=FALSE,lwd=3,bty="n",yaxt="none",font=2, xlim=c(xlim_min,xlim_max),xaxt="none")
    
    mtext(side=1, line=4, paste('Log2 length \n',region,sep=''), col='black', font=2,cex=1.2)
    mtext(side=2, line=3, 'Fn(x)', col="black", font=2, cex=1.2)
    
    axis(side=1,seq(floor(as.numeric(quantile(log2(as.numeric(lenForAnalysis)),0.01))),ceiling(as.numeric(quantile(log2(as.numeric(lenForAnalysis)),0.99))),1), font=2,lwd=2)
    axis(side=2,seq(0,1,0.2), font=2,las=2,lwd=2)
  }
  #plot background
  if(plotType=='violin'){
    vioplot::vioplot(log2(as.numeric(lenForAnalysis)),add=TRUE,at=1,col='grey65', xaxt='n', xlab='', ylab='',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE)
  } else if (plotType=='boxplot'){
    boxplot(log2(as.numeric(lenForAnalysis)),add=TRUE,at=1,col='grey65', xaxt='n', xlab='', ylab='',type='n',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE,outcol='grey65',whiskcol='grey65',outline=FALSE,medcol="black",staplelty = 0,whisklty = 1)
  }
  #
  if(!is.null(geneList)){
    tmpColour <- append(as.character(AnotaColours[regulation]), geneListcolours)
  } else {
    tmpColour <- as.character(AnotaColours[regulation])
  }
  #
  for(j in 1:length(res)){
    tmpVal <- log2(as.numeric(lenForAnalysis[names(lenForAnalysis) %in% res[[j]]]))
    #
    if(plotType=='boxplot'| plotType=='violin'){
      text(j+1,0,round(mean(as.numeric(lenForAnalysis[names(lenForAnalysis) %in% res[[j]]])),0),font=2)
    } else if (plotType=='ecdf'){
      lines(ecdf(tmpVal),col=tmpColour[j],main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=4)
    }
    #
    if(plotType=='violin'){
      vioplot::vioplot(tmpVal,add=TRUE,at=j+1,col=tmpColour[j], xaxt='n', xlab='', ylab='',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE)
    }  else if (plotType=='boxplot'){
      boxplot(tmpVal,add=TRUE,at=j+1,col=tmpColour[j], xaxt='n', xlab='', ylab='',type='n',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE,outcol=tmpColour[j],whiskcol=tmpColour[j],outline=FALSE,medcol="black",staplelty = 0,whisklty = 1)
    }
    resCol[[j+1]] <- tmpVal
  }
  #Plot stats
  if(!is.null(comparisons)){
    if(plotType=='ecdf'){
      tableOut <- matrix(NA, nrow= length(comparisons), ncol= 5)
      colnames(tableOut) <- c('signature','Wilcox_pval','q25','q50','q75')
      colT <- as.character()
    }
    for(j in 1:length(comparisons)){
      compTmp <- comparisons[[j]]+1
      #stats
      pvalTmp <- format(as.numeric(wilcox.test(resCol[[compTmp[1]]], resCol[[compTmp[2]]],alternative='two.sided')[3]),scientific = TRUE,digits=2)
      #
      if(plotType=='boxplot'| plotType=='violin'){
        yposTmp <- range(log2(as.numeric(lenForAnalysis)))[2]+j
        rect(xleft = compTmp[1],xright = compTmp[2],ybottom = yposTmp, ytop = yposTmp,lwd=2)
        #
        text(sum(compTmp)/2,yposTmp+0.5, pvalTmp ,cex=0.75)
      } else if (plotType=='ecdf'){
        #
        colT[j] <-  ifelse(length(which(comparisons[[j]] == 0))>0,as.character(AnotaColours[regulation[comparisons[[j]][comparisons[[j]] != 0]]]),NA)
        #
        if((compTmp[1]-1)==0){
          n1 <- 'bg'
        } else {
          n1 <- paste(paste('c',contrast[(compTmp[1]-1)],sep=''),regulation[(compTmp[1]-1)],sep='_')
        }
        if((compTmp[2]-1)==0){
          n2 <- 'bg'
        } else {
          n2 <- paste(paste('c',contrast[(compTmp[2]-1)],sep=''),regulation[(compTmp[2]-1)],sep='_')
        }
        tableOut[j,1] <- paste(n1,'vs', n2,sep=' ')
        tableOut[j,2] <- pvalTmp
        #Calculate percentiles for Background
        tmpBg <- sort(as.numeric(resCol[[compTmp[1]]]))
        ecdfBg <- 1:length(tmpBg)/length(tmpBg)
        bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
        bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
        bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
        
        #Calculate percentiles for second and difference from background
        tmpSign <- sort(as.numeric(resCol[[compTmp[2]]]))
        ecdfSign <- 1:length(tmpSign)/length(tmpSign)
        tableOut[j,3] <- format(tmpSign[which(ecdfSign >= 0.25)[1]] - bg_025, digits = 2)
        tableOut[j,4] <- format(tmpSign[which(ecdfSign >= 0.5)[1]] - bg_05, digits = 2)
        tableOut[j,5] <- format(tmpSign[which(ecdfSign >= 0.75)[1]] - bg_075, digits = 2)
      }
    }
    if (plotType=='ecdf'){
      plotrix::addtable2plot(xlim_min,1.01,tableOut,bty="n",display.rownames=FALSE,hlines=FALSE,vlines=TRUE,title="",cex = 0.7,bg= colT,xpad=0.1,ypad=1.4,xjust=0,yjust=1)
    }
  }
  dev.off()
  #
  return(lenForAnalysis)
}
