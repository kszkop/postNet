##Run analysis length
contentAnalysis <- function(ads=NULL,
                           regulation=NULL,
                           contrast=NULL,
                           comparisons=NULL,
                           content,
                           region, #UTR5, CDS, UTR3
                           subregion=NULL, #number of nucleotides from start if positive or end if negative.
                           subregionSel, #select or exclude , required if subregion is not null.
                           annot,
                           selection, #shortest, longest, random (default)
                           plotType='boxplot', # option 'violin' or 'ecdf'
                           geneList=NULL,
                           geneListcolours=NULL,
                           customBg=NULL,
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
  if(!is.null(subregion)){
    #
    subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos=subregion,subregionSel=subregionSel)))
  #
    annotBgSel$seqTmp <- subSeq
  }
  annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp),]
  #
  contForAnalysis <- as.numeric()
  for(i in 1:nrow(annotBgSel)){
    tmpSeq <- annotBgSel$seqTmp[i]
    tmpCont <- sapply(seqinr::s2c(toupper(content)), function(x) calc_content(tmpSeq, x))
    contForAnalysis[i] <- sum(tmpCont)
  }
  names(contForAnalysis) <- annotBgSel$geneID
  
  #Prepare plotting
  if(isTRUE(plotOut)){
    #
    resOut <- resSel(vIn=contForAnalysis, ads=ads, regulation=regulation, contrast=contrast, customBg=customBg, geneList=geneList)
    #
    coloursOut <- coloursSel(ads=ads, regulation=regulation, geneList=geneList, geneListcolours=geneListcolours,customBg=customBg)
    #
    #Plot
    pdf(ifelse(is.null(pdfName),paste(region,plotType,'contentAnalysis.pdf',sep='_'), paste(pdfName,region,plotType,'contentAnalysis.pdf',sep='_')),width= 8,height=8, useDingbats = F)
    if(plotType=='boxplot'| plotType=='violin'){
      par(mar=c(8,12,5,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
      #calculate xlim
      if(!is.null(regulation)){
        xlimIn <- c(0.5,length(regulation)+ifelse(!is.null(geneList),length(geneList),0)+1.5)
      } else {
        xlimIn <- c(0.5,length(geneList)+1.5)
      }
      plot(1, 1, xlim=xlimIn, ylim=c(0,100+(length(comparisons)*5)), xaxt='n', xlab='', ylab='',type='n',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE)

      axis(side=2, font=2,las=2,lwd=2,at=seq(0,100,25),labels = seq(0,100,25))
      text(1:length(resOut), par("usr")[3] - 0.45, labels=names(resOut), xpd=NA,cex=0.9,srt=45,adj=1)
      mtext(side=2, line=6, paste(toupper(content), 'content (%)', sep=' '), col="black", font=2, cex=1.7,at=50)
      #
      if(!is.null(ads) | !is.null(customBg)){
        abline(lty=5, h=median(resOut[[1]]))
      }
      for(i in 1:length(resOut)){
        if(plotType=='violin'){
          vioplot::vioplot(resOut[[i]],add=TRUE,at=i,col=coloursOut[i], xaxt='n', xlab='', ylab='',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE)
        } else if (plotType=='boxplot'){
          boxplot(resOut[[i]],add=TRUE,at=i,col=coloursOut[i], xaxt='n', xlab='', ylab='',type='n',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE,outcol='grey65',whiskcol='grey65',outline=FALSE,medcol="black",staplelty = 0,whisklty = 1)
        }
        text(i,0,round(mean(resOut[[i]],0)),font=2)
      }
    } else if (plotType=='ecdf'){
      #
      xlim_min <- as.numeric(quantile(as.numeric(unlist(resOut)),0.01))
      xlim_max <- as.numeric(quantile(as.numeric(unlist(resOut)),0.99))
      #
      par(mar=c(5,5,8,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
      plot(ecdf(resOut[[1]]),col=coloursOut[1],main='',xlab='',ylab='',verticals=TRUE, do.p=FALSE,lwd=3,bty="n",yaxt="none",font=2, xlim=c(xlim_min,xlim_max),xaxt="none")
      
      for(i in 2:length(resOut)){
        lines(ecdf(resOut[[i]]),col=coloursOut[i],main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=4)
      }
      
      mtext(side=1, line=4, paste(paste(toupper(content), 'content (%)', sep=' '),'\n',region,sep=''), col='black', font=2,cex=1.2)
      mtext(side=2, line=3, 'Fn(x)', col="black", font=2, cex=1.2)
      
      axis(side=1,seq(floor(xlim_min),ceiling(xlim_max),10), font=2,lwd=2)
      axis(side=2,seq(0,1,0.2), font=2,las=2,lwd=2)
    }
    #Plot stats
    if(!is.null(comparisons)){
      for(j in 1:length(comparisons)){
        if(!is.null(ads) | !is.null(customBg)){
          compTmp <- comparisons[[j]]+1
        } else {
          compTmp <- comparisons[[j]]
        }
        #stats
        pvalTmp <- format(as.numeric(wilcox.test(resOut[[compTmp[1]]], resOut[[compTmp[2]]],alternative='two.sided')[3]),scientific = TRUE,digits=2)
        #
        if(plotType=='boxplot'| plotType=='violin'){
          yposTmp <- range(as.numeric(unlist(resOut)))[2]+(j*5)
          rect(xleft = compTmp[1],xright = compTmp[2],ybottom = yposTmp, ytop = yposTmp,lwd=2)
          #
          text(sum(compTmp)/2,yposTmp+2.5, pvalTmp ,cex=0.75)
        } else if (plotType=='ecdf'){
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
          plotrix::addtable2plot(xlim_min,1.01,tableOut,bty="n",display.rownames=FALSE,hlines=FALSE,vlines=TRUE,title="",cex = 0.7,bg=colT,xpad=0.1,ypad=1.4,xjust=0,yjust=1)
        }
      }
    }
    dev.off()
  }
  #
  return(contForAnalysis)
}

  
  