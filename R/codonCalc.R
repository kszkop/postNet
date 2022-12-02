#Function based on given AA/codons calculate freq across all genes
codonCalc  <- function(codonIn, #output of codonUsage function
                       analysis, #whether 'AA' or 'codon'
                       featsel, #vector of selected codons
                       ads=NULL,
                       regulation=NULL,
                       contrast=NULL,
                       comparisons=NULL,
                       unit='count', #option 'freq'
                       pdfName=NULL,
                       geneList=NULL,
                       geneListcolours=NULL,
                       customBg=NULL,
                       plotOut=TRUE
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
  
  ##Prepare plotting
  if(isTRUE(plotOut)){
    #
    resOut <- resSel(vIn=codonCalcOut, ads=ads, regulation=regulation, contrast=contrast, customBg=customBg, geneList=geneList)
    #
    coloursOut <- coloursSel(ads=ads, regulation=regulation, geneList=geneList, geneListcolours=geneListcolours,customBg=customBg)
    #
    pdf(ifelse(is.null(pdfName),'codonCalc.pdf', paste(pdfName,'codonCalc.pdf',sep='_')),width= 8,height=8, useDingbats = F)
    par(mar=c(5,5,8,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
    xlim_min <- as.numeric(quantile(as.numeric(unlist(resOut)),0.01))
    xlim_max <- as.numeric(quantile(as.numeric(unlist(resOut)),0.99))
    #
    plot(ecdf(resOut[[1]]),col=coloursOut[1],main='',xlab='',ylab='',verticals=TRUE, do.p=FALSE,lwd=3,bty="n",yaxt="none",font=2, xlim=c(xlim_min,xlim_max),xaxt="none")
    
    for(i in 2:length(resOut)){
      lines(ecdf(resOut[[i]]),col=coloursOut[i],main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=4)
    }
    
    mtext(side=1, line=4, paste('codon usage \n',unit,sep=''), col='black', font=2,cex=1.2)
    mtext(side=2, line=3, 'Fn(x)', col="black", font=2, cex=1.2)
    
    axis(side=1,seq(floor(xlim_min),ceiling(xlim_max),1), font=2,lwd=2)
    axis(side=2,seq(0,1,0.2), font=2,las=2,lwd=2)
    
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
  return(codonCalcOut)
}
