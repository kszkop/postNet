###fold energy
foldEnergyAnalysis <- function(source='load',#option to 'load' available or create new 'create' or provide 'custom'
                                  version,
                                  species,
                                  region=NULL, #UTR5, CDS, UTR3
                                  fromFasta=FALSE,
                                  customFile, #for custom option, path to the file,
                                  onlyRun=FALSE,#if only prepare foldenery file without analysis
                                  ads=NULL,
                                  regulation=NULL,
                                  contrast=NULL,
                                  comparisons=NULL,
                                  geneList=NULL,
                                  geneListcolours=NULL,
                                  annot,
                                  selection, #shortest, longest, random (default)
                                  plotType='boxplot',# option or 'ecdf'
                                  resid=TRUE,
                                  pdfName=NULL,
                                  plotOut=TRUE
){
  if(source=='create'){
    if(isTRUE(fromFasta)){
      #
      runMfold(customFile)
      #
      if(!isTRUE(onlyRun)){
        energyIn <- read.delim(gsub('.fa','_foldEnergy.txt',customFile), stringsAsFactors=FALSE)
        energyIn <- energyIn[!grepl('NM_',energyIn$fold_energy),]
        energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
      }
    } else {
      #Select region of interest
      if(region=='UTR5'){
        #Write out sequences
        seqinr::write.fasta(sequences=as.list(annot$UTR5_seq),names=annot$id,file.out=paste(region,'.fa',sep=''))
        #
        runMfold(paste(region,'.fa',sep=''))
      }
      if(region=='UTR3'){
        #Write out sequences
        seqinr::write.fasta(sequences=as.list(annot$UTR3_seq),names=annot$id,file.out=paste(region,'.fa',sep=''))
        #
        runMfold(paste(region,'.fa',sep=''))
      }
      if(region=='CDS'){
        #Write out sequences
        seqinr::write.fasta(sequences=as.list(annot$CDS_seq),names=annot$id,file.out=paste(region,'.fa',sep=''))
        #
        runMfold(paste(region,'.fa',sep=''))
        #
      }
      #
      if(!isTRUE(onlyRun)){
        energyIn <- read.delim(paste(region,'foldEnergy.txt', sep='_'), stringsAsFactors=FALSE)
        energyIn <- energyIn[!grepl('NM_',energyIn$fold_energy),]
        energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
      }
    }
  } else if(source=='load'){
    #list existing species
    currTmp <- list.files(system.file("extdata/annotation/refseq/",package = "anota2seqUtils"))
    
    if(!species %in% currTmp){
      stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
    }
    if(species=='human'){
      energyIn <- read.delim(system.file(paste("extdata/annotation/refseq/human",version,sep='/'), paste('humanDB_',region,"_foldEnergy", ".txt.gz",sep=''), package = "anota2seqUtils"), stringsAsFactors=FALSE)
      energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
    }
    if(species=='mouse'){
      energyIn <- read.delim(system.file(paste("extdata/annotation/refseq/mouse",version,sep='/'),  paste('mouseDB_',region,"_foldEnergy", ".txt.gz",sep=''), package = "anota2seqUtils"), stringsAsFactors=FALSE)
      energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
    }
  } else if(source=='custom'){
    #
    energyIn <- read.delim(customFile, stringsAsFactors=FALSE)
    energyIn$fold_energy <- as.numeric(energyIn$fold_energy)
  } else {
    stop("No correct option for source file provided")
  }
  #
  if(!isTRUE(onlyRun)){
    #Combine with annot
    energyInGene <- merge(energyIn, annot[,c(1,2)],by='id',all.x = T)
    energyInGene <- na.omit(energyInGene)
    
    #Select per gene level
    if(selection=='shortest'){
      energyInGeneSel <- as.data.frame(energyInGene %>% group_by(geneID) %>% dplyr::slice(which.min(length)))
    } else if(selection=='longest'){
      energyInGeneSel <- as.data.frame(energyInGene %>% group_by(geneID) %>% dplyr::slice(which.max(length)))
    } else {
      energyInGeneSel <- as.data.frame(energyInGene %>% group_by(geneID) %>% dplyr::slice_sample(n = 1))
    }
    #
    if(isTRUE(resid)){
      feForAnalysis <- lm(as.numeric(energyInGeneSel$fold_energy) ~ as.numeric(energyInGeneSel$length))$residuals
      names(feForAnalysis) <- energyInGeneSel$geneID
    } else {
      feForAnalysis <- as.numeric(energyInGeneSel$fold_energy)
      names(feForAnalysis) <- energyInGeneSel$geneID
    }
    #
    if(!is.null(ads)){
      bg <- row.names(ads@dataP)
      if(!is.null(geneList)){
        bg <- unique(c(bg, as.character(unlist(geneList))))
      }
      feForAnalysisBg <- feForAnalysis[names(feForAnalysis) %in% bg]
    } else {
      if(!is.null(customBg)){
        #add here to be sure that all genelist are in bg
        bg <- customBg
        if(!is.null(geneList)){
          bg <- unique(c(bg, as.character(unlist(geneList))))
        }
        feForAnalysisBg <- feForAnalysis[names(feForAnalysis) %in% bg]
      } else {
        feForAnalysisBg <- feForAnalysis[names(feForAnalysis) %in% as.character(unlist(geneList))]
      }
    }
    #
    #Prepare plotting
    if(isTRUE(plotOut)){
      #
      resOut <- resSel(vIn=feForAnalysisBg, ads=ads, regulation=regulation, contrast=contrast, customBg=customBg, geneList=geneList)
      #
      coloursOut <- coloursSel(ads=ads, regulation=regulation, resOut=resOut, geneList=geneList, geneListcolours=geneListcolours,customBg=customBg)
      #
    
      #Plot
      pdf(ifelse(is.null(pdfName),paste(ifelse(!is.null(region),region,''),plotType,'foldEnergyAnalysis.pdf',sep='_'), paste(pdfName,ifelse(!is.null(region),region,''),plotType,'foldenergyAnalysis.pdf',sep='_')),width= 8,height=8, useDingbats = F)
      #
      #c(-roundUpNice(abs(as.numeric(quantile(as.numeric(feForAnalysis),0.01)))), roundUpNice(as.numeric(quantile(as.numeric(feForAnalysis),0.99)))+(length(comparisons)*10))
      #
      xlim_set <- roundUpNice(abs(max(abs(as.numeric(quantile(as.numeric(unlist(resOut)),0.01))),abs(as.numeric(quantile(as.numeric(unlist(resOut)),0.99))))))
      #
      if(plotType=='boxplot'|plotType=='violin'){
        #calculate xlim
        if(!is.null(regulation)){
          xlimIn <- c(0.5,length(regulation)+ifelse(!is.null(geneList),length(geneList),0)+1.5)
        } else {
          xlimIn <- c(0.5,length(geneList)+1.5)
        }
        par(mar=c(8,5,8,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
        plot(1, 1, xlim=xlimIn, ylim=c(-xlim_set, ifelse(isTRUE(resid),xlim_set,50)), xaxt='n', xlab='', ylab='',type='n',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE)
        for(i in 1:length(resOut)){
          if(plotType=='violin'){
            vioplot::vioplot(resOut[[i]],add=TRUE,at=i,col=coloursOut[i], xaxt='n', xlab='', ylab='',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE)
          } else if (plotType=='boxplot'){
            boxplot(resOut[[i]],add=TRUE,at=i,col=coloursOut[i], xaxt='n', xlab='', ylab='',type='n',main="",lwd=1,bty="n",yaxt="n",font=2,frame.plot=FALSE,outcol='grey65',whiskcol='grey65',outline=FALSE,medcol="black",staplelty = 0,whisklty = 1)
          }
          text(i,-xlim_set,round(mean(resOut[[i]],0)),font=2)
        }
        axis(side=2, font=2,las=2,lwd=2,at=seq(-xlim_set,ifelse(isTRUE(resid),xlim_set,50),20),labels = seq(-xlim_set,ifelse(isTRUE(resid),xlim_set,50),20))
      
        mtext(side=2, line=6, 'fold energy', col="black", font=2, cex=1.7,at=0)
        if(!is.null(ads) | !is.null(customBg)){
          abline(lty=5, h=median(resOut[[1]]))
        }
        text(1:length(resOut), par("usr")[3] - 0.45, labels=names(resOut), xpd=NA,cex=0.9,srt=45,adj=1)
        
      } else if (plotType=='ecdf'){
        par(mar=c(5,5,8,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.4, cex.main=1.7,cex.lab=1.3)
        plot(ecdf(resOut[[1]]),col=coloursOut[1],main='',xlab='',ylab='',verticals=TRUE, do.p=FALSE,lwd=3,bty="n",yaxt="none",font=2, xlim=c(-xlim_set, ifelse(isTRUE(resid),xlim_set,0)) ,xaxt="none")
      
        mtext(side=1, line=4, paste('fold energy','\n',paste(ifelse(!is.null(region),region,''),sep='')), col='black', font=2,cex=1.2)
        mtext(side=2, line=3, 'Fn(x)', col="black", font=2, cex=1.2)
      
        axis(side=1,seq(-xlim_set,ifelse(isTRUE(resid),xlim_set,0),20), font=2,lwd=2)
        axis(side=2,seq(0,1,0.2), font=2,las=2,lwd=2)
        
        for(i in 2:length(resOut)){
          lines(ecdf(resOut[[i]]),col=coloursOut[i],main='',xlab='',verticals=TRUE, do.p=FALSE,lwd=4)
        }
      }
      #
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
    return(feForAnalysisBg)
  }
}
