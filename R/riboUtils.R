riboUtils <- function(annot, adjustStart = 0, fastaFile, path=NULL, name=NULL, runRUST=TRUE){
  #
  if(!is.null(path)){
    if(unlist(strsplit(path,''))[length(unlist(strsplit(path,'')))]!='/'){
      path <- paste(path,'/',sep='')      
    }
  }
  #
  seq <- list()
  seq[[1]] <- annot$UTR5_seq
  lenUTR5 <- as.numeric(sapply(annot$UTR5_seq, function(x) length(seqinr::s2c(x))))
  seq[[2]] <- annot$CDS_seq
  lenCDS <- as.numeric(sapply(annot$CDS_seq, function(x) length(seqinr::s2c(x))))
  seq[[3]] <- annot$UTR3_seq
  lenUTR3 <- as.numeric(sapply(annot$UTR3_seq, function(x) length(seqinr::s2c(x))))
  
  #
  seqAll <- combSeq(seqIn = seq)
  lenTmp <- as.numeric(sapply(seqAll, function(x) length(seqinr::s2c(x))))
  #
  annot$lenTmp <- lenTmp
  annot$lenUTR5 <- lenUTR5
  annot$lenCDS <- lenCDS
  annot$lenUTR3 <- lenUTR3
  #
  annotSel <- isoSel(annot=annot, method='longest')
  
  #ribowaltz
  annotRW <- annot[,c(1,6:9)]
  colnames(annotRW) <- c('transcript','l_tr','l_utr5','l_cds','l_utr3')
  annotRW <- data.table::data.table(annotRW)
  
  #ID:
  id <- list.files(ifelse(is.null(path),'.',path), pattern='.bam$')
  reads_list <- riboWaltz::bamtolist(bamfolder = ifelse(is.null(path),'.',path), annotation = annotRW)

  # Filter the reads based on a certain length
  rLen <- lapply(reads_list,function(x) pull(x,var=4))
  rLen <- unlist(rLen,recursive = FALSE, use.names = FALSE)
  rLen <- table(rLen)
  rSel <- names(which(rLen/sum(rLen)>0.075))
  #Remove everything below 25 or above 35 as these are probably not really riboseq reads anyway
  rSel <- rSel[which(rSel>24 & rSel < 36)]
  #
  filtered_list <- riboWaltz::length_filter(data = reads_list, length_filter_mode = "custom",length_range = min(rSel):max(rSel))
  # Caluclate P-site information based on the filtered lengths
  pInfoOut <- riboWaltz::psite(filtered_list,plot=T,plot_dir = ".",plot_format = "pdf")
  reads_pSite_info <- riboWaltz::psite_info(filtered_list, offset = pInfoOut)
  
  #Run thorough samples
  countsOut <- list()
  for(i in 1:length(reads_pSite_info)){
    #Extract data for the sample
    idFile <- id[i]
    #extract data for that file
    sTmp <- reads_pSite_info[[i]]
    
    #Transcript level
    sTmp_trans <- sTmp[,c(1,3)]
    #Count
    counts_trans <- sTmp_trans %>% group_by(transcript, psite) %>% summarise(count=n())
    #Convert to dataList format
    dataList_trans <- lapply(split(counts_trans[-1], as.character(counts_trans[[1]])), function(x) setNames(x$count, x$psite))
    cat("\n","Writing out transcript level dpn","\n")
    outPut_trans <- as.matrix(unlist(sapply(unique(names(dataList_trans)), generateOut, dataList_trans)))
    write.table(outPut_trans, ifelse(is.null(path),paste(gsub('.bam','',idFile), '_trans.dpn',sep=''),paste(path, paste(gsub('.bam','',idFile), 'trans.dpn',sep='_'),sep='/')), sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
    
    #Gene level
    #first subset reads only to the selected isoform.
    sTmp_sel <- sTmp[transcript %in% annotSel$id]
    #adjust position if desired (maybe better to adjust position of cds and recalculate 
    if(adjustStart != 0){
      #
      cdsStartTmp <- sTmp_sel$cds_start
      cdsStartAdj <- cdsStartTmp + adjustStart
      #
      sTmp_sel$cds_start <- cdsStartAdj
      #but make sure the start is not after the end after adjustment as it might causes stragne things). I guess I will remove these sitations
      cdsStopTmp <- sTmp_sel$cds_stop
      checkCDSadj <- cdsStartAdj - cdsStopTmp
      if(length(which(checkCDSadj >= 0))>0){
        transTorem <- as.character(sTmp_sel$transcript)[which(checkCDSadj >= 0)]
        sTmp_sel <- sTmp_sel[!as.character(sTmp_sel$transcript) %in% transTorem,]
      }
      #recalculate region
      psite_from_startTmp <- sTmp_sel$psite - sTmp_sel$cds_start
      psite_from_stopTmp <- sTmp_sel$psite - sTmp_sel$cds_stop
      #
      psite_regionTmp <- rep('NA', length(psite_from_startTmp ))
      #
      psite_regionTmp[which(psite_from_startTmp >= 0 & psite_from_stopTmp <= 0 )] <- 'cds'
      psite_regionTmp[which(psite_from_startTmp < 0)] <- '5utr'
      psite_regionTmp[which(psite_from_stopTmp > 0)] <- '3utr'
      #
      sTmp_sel$psite_from_start <- psite_from_startTmp
      sTmp_sel$psite_from_stop <- psite_from_stopTmp
      sTmp_sel$psite_region <- psite_regionTmp
    }
    #add gene IDs
    geneInfo <- annotSel[,c(1,2)]
    colnames(geneInfo)[1] <- 'transcript'
    sTmp_sel[geneInfo, on = 'transcript', geneID := i.geneID]
    sTmp_sel <- sTmp_sel[sTmp_sel$psite_region=='cds']
    
    #Count
    counts_geneTmp <- sTmp_sel$geneID
    counts_gene <- table(counts_geneTmp)
    #
    countsOut[[gsub('.bam','',idFile)]] <- counts_gene
    
    #Ribowatlz
    pdf(gsub('.bam','_riboWaltzQC_Filt.pdf',idFile),height=18,width=12)
    gridExtra::grid.arrange(riboWaltz::rlength_distr(reads_list,sample=gsub('.bam','',idFile),cl=99)[[paste("plot",gsub('.bam','',idFile),sep='_')]],
                            riboWaltz::rends_heat(filtered_list,annotRW,sample=gsub('.bam','',idFile),cl=85,utr5l = 25, cdsl = 40, utr3l = 25)[["plot"]],
                            riboWaltz::region_psite(reads_pSite_info,annotRW,sample=gsub('.bam','',idFile))[["plot"]],
                            riboWaltz::frame_psite_length(reads_pSite_info, sample = gsub('.bam','',idFile),region = "all", cl = 90)[["plot"]],
                            riboWaltz::frame_psite(reads_pSite_info, sample = gsub('.bam','',idFile), region = "all")[["plot"]],
                            riboWaltz::codon_usage_psite(reads_pSite_info, annotRW, sample = gsub('.bam','',idFile),fastapath = fastaFile,fasta_genome = FALSE,frequency_normalization = FALSE)[["plot"]],
                 ncol=2,nrow=3,widths=c(0.33,0.67)
    )
    dev.off()
  
    if(isTRUE(runRUST)){
      ###RUST - maybe convert to R one day if authors agree.
      # Extract dominant read length
      rLenTmp <- as.numeric(names(which.max(table(sTmp$length))))
      offsetLen <- sTmp_sel[sTmp_sel$length==rLenTmp,]$psite[1]-sTmp_sel[sTmp_sel$length==rLenTmp,]$end5[1]
      #
      nameTmpIn1 <- ifelse(is.null(path),idFile, paste(path,idFile,sep='/'))
      nameTmpOut1 <-ifelse(is.null(path),paste(gsub('.bam','',idFile), "_length",rLenTmp,".bam",sep=''),paste(path,paste(gsub('.bam','',idFile), "_length",rLenTmp,".bam",sep=''),sep='/'))
      command1 <- paste("samtools view -h -F 0x10", nameTmpIn1, "| awk '(length($10) ==",rLenTmp,") || $1 ~ /^@/' | samtools view -bS - >", nameTmpOut1, sep=" ")
      system(command1)
      #
    
      command2 <- paste("samtools index",nameTmpOut1, sep=" ")
      system(command2)
      #
      command3 <- paste("python3 rust_codon_p3.py", fastaFile, nameTmpOut1, offsetLen, paste(min(rSel),max(rSel),sep=':'),sep=' ')
      system(command3)
    
      ####
      profile <- read.csv(paste("codon/RUST_codon_file_",paste(gsub('.bam','',idFile), "_length",rLenTmp,".bam_",paste(offsetLen,min(rSel),max(rSel),sep='_'),sep=''),sep=''))
      #
      codPlot(rust=profile,name=paste('codon/RUSTcodon_profile',idFile,sep='_'))
    }
  }
  #combine counts
  countsGene <- t(plyr::ldply(countsOut , rbind, .id=NULL))
  colnames(countsGene) <- names(countsOut)
  #Replace NA with 0
  countsGene[is.na(countsGene)]<-0
  #
  countsGene <- data.frame(geneID=row.names(countsGene), countsGene, row.names = NULL)
  #
  write.table(countsGene,file=ifelse(is.null(name),'riboSeq_counts_genelevel.txt',paste(name,'riboSeq_counts_genelevel.txt',sep='_')), col.names=T,row.names=F,sep='\t',quote=F)
}

