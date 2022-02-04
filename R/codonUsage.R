##Run analysis length
codonUsage <- function(ads,
                          analysis,
                          regulation,
                          contrast,
                          type='sequence',#option: 'sequence' or if ribosomal profiling 'riboSeq'
                          comparisons=NULL,
                          annotType='ccds',#option: 'refseq' or 'ccds'
                          annot=NULL, #it is required for refseq, if ccds given leave it null. and give options for sourse, species,
                          source='load',#option to 'load' available or create new 'create'or provide 'custom'
                          species=NULL, #it is required for ccds
                          subregion=NULL, #number of nucleotides from start if positive or end if negative.
                          subregionSel, #select or exclude , required if subregion is not null.
                          geneList=NULL,
                          geneListnames=NULL,
                          geneListcolours=NULL,
                          selection, #shortest, longest, random (default)
                          pAdj=0.01,
                          plotHeatmap=TRUE,
                          pdfName=NULL
){
  #
  if(annotType=='ccds'){
    if(source=='create'){
      if(species=='human'){
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt", destfile = "CCDS_human.txt")
        ccds <- read.delim(file = "CCDS_human.txt", as.is = c(F, T, T, F, T, F, F, T, T, T, F))
        unlink("CCDS_human.txt")
        #
        ccds$ccds_id_cl <-  gsub("\\..*","",ccds$ccds_id )
        
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS_nucleotide.current.fna.gz", 
                      destfile = "CCDS_nucleotide_human.fna.gz")
        gunzip("CCDS_nucleotide_human.fna.gz")
        ccdsSeq <- seqinr::read.fasta("CCDS_nucleotide_human.fna",seqtype = 'AA')
        unlink("CCDS_nucleotide_human.fna")
        
        # Below, I just exclude these duplicated sequences 
        ccdsSeq <- ccdsSeq[!duplicated(sapply(strsplit(names(ccdsSeq), "\\|"), function(x) x[1]))]
        
        ccdsSeq <- data.frame(id=sub("\\..*", "", names(ccdsSeq)), CDS_seq=t(as.data.frame(lapply(ccdsSeq,function(x) paste(x, collapse='')))), row.names=NULL, stringsAsFactors = FALSE)
        
        annot <- merge(unique(ccds[, c("ccds_id_cl", "gene")]), ccdsSeq, by.x = "ccds_id_cl", by.y = "id", all.y = TRUE, all.x = FALSE)
        colnames(annot) <- c('id','geneID','CDS_seq')
        write.table(annot, file='customDB_ccds.txt',col.names=T, row.names=F, sep='\t',quote=F)
        
      } else if(species=='mouse'){
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_mouse/CCDS.current.txt", destfile = "CCDS_mouse.txt")
        ccds <- read.delim(file = "CCDS_mouse.txt", as.is = c(F, T, T, F, T, F, F, T, T, T, F))
        unlink("CCDS_mouse.txt")
        #
        ccds$ccds_id_cl <-  gsub("\\..*","",ccds$ccds_id )
        
        #
        download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_mouse/CCDS_nucleotide.current.fna.gz", 
                      destfile = "CCDS_nucleotide_mouse.fna.gz")
        gunzip("CCDS_nucleotide_mouse.fna.gz")
        ccdsSeq <- seqinr::read.fasta("CCDS_nucleotide_mouse.fna",seqtype = 'AA')
        unlink("CCDS_nucleotide_mouse.fna")
        
        # Below, I just exclude these duplicated sequences 
        ccdsSeq <- ccdsSeq[!duplicated(sapply(strsplit(names(ccdsSeq), "\\|"), function(x) x[1]))]
        
        ccdsSeq <- data.frame(id=sub("\\..*", "", names(ccdsSeq)), CDS_seq=t(as.data.frame(lapply(ccdsSeq,function(x) paste(x, collapse='')))), row.names=NULL, stringsAsFactors = FALSE)
        
        annot <- merge(unique(ccds[, c("ccds_id_cl", "gene")]), ccdsSeq, by.x = "ccds_id_cl", by.y = "id", all.y = TRUE, all.x = FALSE)
        colnames(annot) <- c('id','geneID','CDS_seq')
        write.table(annot, file='customDB_ccds.txt',col.names=T, row.names=F, sep='\t',quote=F)
      }
    } else if(source=="load") {
      #list existing species
      currTmp <- list.files(system.file("extdata/annotation/ccds",package = "anota2seqUtils"))
      
      if(!species %in% currTmp){
        stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
      }
      if(species=='human'){
        annot <- read.delim(system.file(paste("extdata/annotation/ccds/human",sep='/'), "humanDB_ccds.txt.gz", package = "anota2seqUtils"), stringsAsFactors=FALSE)
      }
      if(species=='mouse'){
        annot <- read.delim(system.file(paste("extdata/annotation/ccds/mouse",sep='/'), "mouseDB_ccds.txt.gz", package = "anota2seqUtils"), stringsAsFactors=FALSE)  # }
      }
    } else if (source=="custom"){
      annot <- read.delim(customFile, stringsAsFactors=FALSE)
    } else {
      stop("No correct option for annotation file provided")
    }
  }
  #Subset annot for only expressed genes
  bg <- row.names(ads@dataP)
  annotBg <- annot[annot$geneID %in% bg,]
  
  # 
  seqTmp <- annotBg$CDS_seq
  lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
  annotBg <- cbind(annotBg[,c(1:2)], seqTmp, lenTmp)

  #Select per gene level
  if(selection=='shortest'){
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice(which.min(lenTmp)))
  } else if(selection=='longest'){
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice(which.max(lenTmp)))
  } else {
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice_sample(n = 1))
  }
  
  if(!is.null(subregion)){
    #
    subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos=subregion,subregionSel=subregionSel)))
    #
    annotBgSel$seqTmp <- subSeq
  }
  annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp),]
  
  if(type=='sequence'){
    codonTmp <- list()
    for(i in 1:nrow(annotBgSel)){
      codonTmp[[i]] <- codonCount(gene=annotBgSel$geneID[i], seq=annotBgSel$seqTmp[i])
    }
    codonAll <- do.call(rbind, lapply(codonTmp, data.frame, stringsAsFactors=FALSE))
    codonAll <- codonAll[!codonAll$AA %in% c("Stp", "Trp", "Met"),]
  } else if(type=='riboSeq'){
      stop("This option is not yet fully functional")
  }
  #gene list instead of ads object
  if(!is.null(geneList)){
    res <- geneList
    names(res) <- geneListnames
  } else {
    #Extract all results
    results <- anota2seqGetDirectedRegulations(ads)
  
    res <- vector("list", length = length(regulation))
    for(i in unique(contrast)){
      resTmp <- results[[i]][names(results[[i]]) %in% regulation[which(contrast==i)]]
      res[which(contrast==i)] <- resTmp
    }
  }
  if(length(res) < 2){
    stop(
      "For comparison of codon composition between regulations, at least two regulations should be provided.",
      call. = FALSE
    )
  }
  
  #
  resOut <- list()
  for(j in 1:length(res)){
    tmp <- codonAll[codonAll$geneID %in% res[[j]],]
    #
    if(analysis=='codon'){
      tmp <- tmp %>% group_by(codon) %>% summarise(codonPerReg=sum(codonCount))
      #
      tmpSum <- tmp$codonPerReg
      names(tmpSum) <- tmp$codon
    } else if(analysis=='AA'){
      #
      tmp <- tmp %>% group_by(AA) %>% summarise(AAPerReg=sum(codonCount))
      #
      tmpSum <- tmp$AAPerReg
      names(tmpSum) <- tmp$AA
    }
    #
    if(!is.null(geneList)){
      resOut[[geneListnames[j]]] <-  tmpSum
    } else {
      resOut[[regulation[j]]] <-  tmpSum
    }
  }
  resOut <- do.call(rbind, resOut)
  
  #
  chisqTest <- chisq.test(resOut)
  
  if(is.na(chisqTest$p.value)){
    stop("The Chi-squared test could not be performed. Please check that all codons have at least one count for at least one geneset.")
  }  
      
  if(sum(apply(resOut, 2, min) < 5) > 0){
    warning("In the contingency table (counts of codons by geneset), some counts are lower than 5 which may invalidate the chi-squared test. It might be relevant to perform the analysis on a subset or groups of codons instead",
            call. = FALSE)
  }
  
  StResiduals <- chisqTest$stdres
  signifLimit <- sqrt(qchisq(pAdj/(dim(resOut)[1]*dim(resOut)[2]), lower.tail = FALSE, df = 1))
  
  if(isTRUE(plotHeatmap)){
    #
    residOut <- as.matrix(StResiduals)
    if(!is.null(geneList)){
      rownames(residOut) <- geneListnames
    } else {
      rownames(residOut) <- regulation
    }
    residOut  <- t(residOut)
    
    pdf(ifelse(is.null(pdfName),paste(analysis,'Heatmap.pdf',sep=''), paste(pdfName,paste(analysis,'Heatmap.pdf',sep=''),sep='_')),width= 20,height=28, useDingbats = F)
    par(mar=c(10,5,5,5),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.9, cex.main=0.7,cex.lab=0.7)
    col <- grDevices::colorRampPalette(c("blue","white", "red"))(50)
    gplots::heatmap.2(residOut,
              col          = col,
              breaks       = c(seq(-ceiling(max(abs(residOut[,1]))),ceiling(max(abs(residOut[,1]))),length=51)),
              margins      = c(30, 50),
              key =TRUE,
              keysize      = 0.5,
              dendrogram   = "both", 
              trace        = "none", 
              density.info = "none",
              labCol       = colnames(residOut),
              key.par = list(cex=0.90),
              cexCol       = 2.1,
              cexRow       = 1.4,
              key.xlab     = "",
              lhei=c(3,25), 
              lwid=c(3,13),
              na.rm=TRUE,
              main=''
    )
    dev.off()
  }
  #
  signSel <- colnames(StResiduals)[which(abs(StResiduals[1,]) > signifLimit)]
  #
  if(is.null(comparisons)){
    stop(
      "For further steps pairs of comparison of codon composition between regulations should be provided with direction for each.",
      call. = FALSE
    )
  }
  #
  #Prepare comparisons per regulation
  for(j in 1:length(comparisons)){
    compTmp <- comparisons[[j]]
    
    if(analysis=='codon'){
      resTmp <- res[compTmp]
      if(!is.null(geneList)){
        regTmp <- geneListnames[compTmp]
      } else {
        regTmp <- regulation[compTmp]
      }
      #
      compOut <- list()
      for(i in 1:2){
        tmp <- codonAll[codonAll$geneID %in% resTmp[[i]],]
        #  
        tmp <- tmp %>% group_by(codon) %>% summarise(freqPerReg=mean(codonFreq))
        #
        tmpFreq <- tmp$freqPerReg
        names(tmpFreq) <- tmp$codon
      
        #
        compOut[[regTmp[i]]] <-  tmpFreq
      }
      compOut <- as.data.frame(do.call(cbind, compOut))
      compOut$codon <- row.names(compOut)
      compOut$AA <-seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(compOut)))))
      compOut$AAcodon <- paste(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(compOut)))),rownames(compOut),sep='_')

      #Plot frequency
      #
      xylim <- roundUpNice(max(compOut[,c(1,2)])) 
      #
      pdf(ifelse(is.null(pdfName),paste(paste(colnames(compOut)[1:2], collapse='_'),'averageFreq.pdf',sep='_'), paste(pdfName,paste(colnames(compOut)[1:2], collapse='_'),'averageFreq.pdf',sep='_')),width= 8,height=8, useDingbats = F)
      par(mar=c(10,5,5,10),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.9, cex.main=0.7,cex.lab=0.7)
      pOut<- ggplot2::ggplot(compOut,  ggplot2::aes(!!sym(colnames(compOut)[1]),!!sym(colnames(compOut)[2]), col = AA)) +
        ggplot2::theme_bw() +
        ggplot2::xlim(0,xylim) +
        ggplot2::ylim(0,xylim) +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::coord_fixed() +
        ggplot2::geom_line(ggplot2::aes(group = AA), col = "gray", linetype = 1, size = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggrepel::geom_text_repel(ggplot2::aes(label = AAcodon), size = 3, segment.size = 0.2) + 
        ggplot2::labs(x = paste(colnames(compOut)[1],"\n(average codon frequency - per thousand)","\n is it really per thousand ??", sep=''), 
           y = paste(colnames(compOut)[2],"\n(average codon frequency - per thousand)","\n is it really per thousand ??",sep=''))
      print(pOut)
      dev.off()
    
      #
      compOut <- list()
      for(i in 1:2){
        tmp <- codonAll[codonAll$geneID %in% resTmp[[i]],]
        #
        tmp <- tmp %>% group_by(codon) %>% mutate(codonPerReg=sum(codonCount))
        tmp <- tmp %>% group_by(AA) %>% mutate(AAPerReg=sum(AACountPerGene))
      
        #
        tmp <- subset(tmp,!duplicated(codon))
        tmp$codonNormAA <- tmp$codonPerReg/tmp$AAPerReg
      
        #
        tmpcodonNormAA <- tmp$codonNormAA
        names(tmpcodonNormAA) <- tmp$codon
      
        #
        compOut[[regTmp[i]]] <-  tmpcodonNormAA
      }
      compOut <- as.data.frame(do.call(cbind, compOut))
      compOut$codon <- row.names(compOut)
      compOut$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(compOut)))))
      compOut$AAcodon <- paste(seqinr::translate(seqinr::s2c(seqinr::c2s(rownames(compOut)))),rownames(compOut),sep='_')
    
      #Plot codonusage
      #
      xylim <- roundUpNice(max(compOut[,c(1,2)])) 
      #
      pdf(ifelse(is.null(pdfName),paste(paste(colnames(compOut)[1:2], collapse='_'),'codonUsage.pdf',sep='_'), paste(pdfName,paste(colnames(compOut)[1:2], collapse='_'),'codonUsage.pdf',sep='_')),width= 8,height=8, useDingbats = F)
      par(mar=c(10,5,5,10),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.9, cex.main=0.7,cex.lab=0.7)
      pOut <- ggplot2::ggplot(compOut,  ggplot2::aes(!!sym(colnames(compOut)[1]),!!sym(colnames(compOut)[2]), col = AA)) +
        ggplot2::theme_bw() +
        ggplot2::xlim(0,xylim) +
        ggplot2::ylim(0,xylim) +
        ggplot2::geom_point(size = 0.5) +
        ggplot2::coord_fixed() +
        ggplot2::geom_line(ggplot2::aes(group = AA), col = "gray", linetype = 1, size = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggrepel::geom_text_repel(ggplot2::aes(label = AAcodon), size = 3, segment.size = 0.2) + 
        ggplot2::labs(x = paste(colnames(compOut)[1],"\n(Amino acid normalised codon usage)", sep=''), 
             y = paste(colnames(compOut)[2],"\n(Amino acid normalised codon usage)",sep=''))
      print(pOut)
      dev.off()
    #
      resOutTmp <- data.frame(codon= colnames(resOut) , t(resOut[compTmp,]),row.names = NULL)
      #
      resOutTmp$AA <- seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(resOutTmp$codon))))
      #
      if(!is.null(geneList)){
        regs <- geneListnames[compTmp]
      } else {
        regs <- regulation[compTmp]
      }
      regComb <- combn(x=regs,m=2)
    
      statOut <- statOnDf(df = resOutTmp,regs=regComb, analysis = analysis)
      
      #calcuate frequencies per codon
      freqTotal <- list()
      #select genes for both compared categories
      gTmp <- union(resTmp[[1]],resTmp[[2]])
      #
      tmp <- codonAll[codonAll$geneID %in% gTmp,]
      #
      tmp <- tmp %>% group_by(codon) %>% summarise(freq=sum(codonFreq))
      
      #
      finalOut <- merge(tmp, data.frame(statOut), by.x='codon',by.y='row.names')
      #remove not significant from chipseq
      finalOut <- finalOut[finalOut$codon %in% signSel,]
      #
      pdf(ifelse(is.null(pdfName),paste(analysis,'codon_oddratio_freq.pdf',sep=''), paste(pdfName,paste(analysis,'codon_oddratio_freq.pdf',sep=''),sep='_')),width= 8,height=8, useDingbats = F)
      par(mar=c(5,5,5,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.7, cex.main=1.7,cex.lab=1.3)
      plot(finalOut$statOut,finalOut$freq,col='black',pch=20,cex=0.1,xlab='',ylab='',lwd=1,bty="n",xaxt="n",yaxt="n",font=2,xlim=range(finalOut$statOut),ylim=range(finalOut$freq))
      text(finalOut$statOut,finalOut$freq, finalOut$codon,col='black',font=2)
      
      #tmp KS
      #text(finalOut[finalOut$codon %in% c('GAG','CTG','GTG','GAC','GCC'),]$statOut,finalOut[finalOut$codon %in% c('GAG','CTG','GTG','GAC','GCC'),]$freq, finalOut[finalOut$codon %in% c('GAG','CTG','GTG','GAC','GCC'),]$codon,col='brown1',font=2)
      #text(finalOut[finalOut$codon %in% c('AAA','GAA','GAT','AAT','TTT','ATT'),]$statOut,finalOut[finalOut$codon %in% c('AAA','GAA','GAT','AAT','TTT','ATT'),]$freq, finalOut[finalOut$codon %in% c('AAA','GAA','GAT','AAT','TTT','ATT'),]$codon,col='dodgerblue1',font=2)
      #legend(1.5,40, fill=c('brown1','dodgerblue1'),bty='n', c('codonsDown1','codonsDown2'))
      #text(finalOut[finalOut$codon %in% c('AAG','ATC'),]$statOut,finalOut[finalOut$codon %in% c('AAG','ATC'),]$freq, finalOut[finalOut$codon %in% c('AAG','ATC'),]$codon,col='brown1',font=2)
      #legend(0.8,80, fill=c('brown1'),bty='n', c('codonsUp'))
      
      #text(finalOut[finalOut$codon %in% c('CTG','GAG','AAG','CAG'),]$statOut,finalOut[finalOut$codon %in% c('CTG','GAG','AAG','CAG'),]$freq, finalOut[finalOut$codon %in% c('CTG','GAG','AAG','CAG'),]$codon,col='brown1',font=2)
      #text(finalOut[finalOut$codon %in% c('AAA','GAA','GAT','CCT'),]$statOut,finalOut[finalOut$codon %in% c('AAA','GAA','GAT','CCT'),]$freq, finalOut[finalOut$codon %in% c('AAA','GAA','GAT','CCT'),]$codon,col='dodgerblue1',font=2)
      #legend(1.5,40, fill=c('brown1','dodgerblue1'),bty='n', c('codonsUp','codonsDown'))
      
      text(finalOut[finalOut$codon %in% c('CTG','GAG','CAG'),]$statOut,finalOut[finalOut$codon %in% c('CTG','GAG','CAG'),]$freq, finalOut[finalOut$codon %in% c('CTG','GAG','CAG'),]$codon,col='dodgerblue1',font=2)
      text(finalOut[finalOut$codon %in% c('AAA','GAA'),]$statOut,finalOut[finalOut$codon %in% c('AAA','GAA'),]$freq, finalOut[finalOut$codon %in% c('AAA','GAA'),]$codon,col='brown1',font=2)
      legend(0.5,20, fill=c('dodgerblue1','brown1'),bty='n', c('codonsDown1','codonsDown2'))
      
      #tmp KS
      
      mtext(side=2, line=3, 'frequency ', col="black", font=2, cex=1.7)
      mtext(side=1, line=3, 'odd ratio', col="black", font=2, cex=1.7,at=25)
      
      axis(side=2,seq(floor(range(finalOut$freq)[1]),ceiling(range(finalOut$freq)[2]),10), font=2,las=2,lwd=2)
      axis(side=1,seq(floor(range(finalOut$statOut)[1]),ceiling(range(finalOut$statOut)[2]),0.25), font=2,lwd=2)
      dev.off()
    } else if(analysis=='AA'){
      resOutTmp <- data.frame(AA= colnames(resOut) , t(resOut[compTmp,]),row.names = NULL)
      
      #
      if(!is.null(geneList)){
        regs <- geneListnames[compTmp]
      } else {
        regs <- regulation[compTmp]
      }
      regComb <- combn(x=regs,m=2)
      
      statOut <- statOnDf(df = resOutTmp,regs=regComb, analysis = analysis)

      #calcuate frequencies per AA
      freqTotal <- list()
      #select genes for both compared categories
      gTmp <- union(resTmp[[1]],resTmp[[2]])
      #
      tmp <- codonAll[codonAll$geneID %in% gTmp,]
      #
      tmp <- tmp %>% group_by(AA) %>% summarise(freq=sum(codonFreq))
      
      #
      finalOut <- merge(tmp, data.frame(statOut), by.x='AA',by.y='row.names')
      #remove not significant from chipseq
      finalOut <- finalOut[finalOut$AA %in% signSel,]
      #
      
      pdf('AA_oddratio_freq.pdf',width= 8,height=8, useDingbats = F)
      par(mar=c(5,5,5,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=1.7, cex.main=1.7,cex.lab=1.3)
      plot(finalOut$statOut,finalOut$freq,col='black',pch=20,cex=0.1,xlab='',ylab='',lwd=1,bty="n",xaxt="n",yaxt="n",font=2,xlim=range(finalOut$statOut),ylim=range(finalOut$freq))
      text(finalOut$statOut,finalOut$freq, finalOut$AA,col='black')
      
      mtext(side=2, line=3, 'frequency ', col="black", font=2, cex=1.7)
      mtext(side=1, line=3, 'odd ratio', col="black", font=2, cex=1.7,at=25)
      
      axis(side=2,seq(floor(range(finalOut$freq)[1]),ceiling(range(finalOut$freq)[2]),10), font=2,las=2,lwd=2)
      axis(side=1,seq(floor(range(finalOut$statOut)[1]),ceiling(range(finalOut$statOut)[2]),0.25), font=2,lwd=2)
      dev.off()
    }
  }
  return(codonAll)
}


  
  
 
  
  
  
  
  
  
  