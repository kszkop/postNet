#Run motif analysis
motifAnalysis <- function(ads,
                          annot,
                          regulation=NULL,
                          contrast=NULL,
                          geneVec=NULL,
                          geneVecName=NULL,
                          region, #UTR5, CDS, UTR3
                          subregion=NULL, #number of nucleotides from start if positive or end if negative.
                          subregionSel=NULL, #select or exclude , required if subregion is not null.
                          selection, #shortest, longest, random (default)
                          stremeThreshold = 0.05, #pvalue , max is 0.05 
                          minwidth=6, #min motif width (default for original software is 8 but thought I would reduce it a bit as 8 sounds bit harsh)
                          stremeName=NULL, #name for output folder
                          tomtomName=NULL,
                          tomtom_database = NULL#path to .meme database of the choice. dowload from website
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
  } else if(selection=='longest'){
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice(which.max(lenTmp)))
  } else {
    annotBgSel <- as.data.frame(annotBg %>% group_by(geneID) %>% dplyr::slice_sample(n = 1))
  }
  #
  if(!is.null(subregion)){
    #
    subSeq <- as.character(sapply(annotBgSel$seqTmp, function(x) subset_seq(x, pos=subregion,subregionSel=subregionSel)))
    #
    annotBgSel$seqTmp <- subSeq
  }
  annotBgSel <- annotBgSel[!is.na(annotBgSel$seqTmp),]
  #
  seqForAnalysis <- annotBgSel$seqTmp
  names(seqForAnalysis) <- annotBgSel$geneID

  #Select regualted genes
  if(is.null(geneVec)){
    #Extract results
    results <- anota2seqGetDirectedRegulations(ads)
    geneSel <- results[[contrast]][[regulation]]
  } else {
    geneSel <- geneVec
  }
  #create vector for regulated and control sequences
  regSeq <- seqForAnalysis[names(seqForAnalysis) %in% geneSel]
  controlSeq <- seqForAnalysis#[!names(seqForAnalysis) %in% geneVec]
  
  #write fasta
  seqinr::write.fasta(sequences=as.list(as.character(regSeq)),names=names(regSeq),file.out=ifelse(is.null(geneVec),paste(paste('Regulated_contrast',contrast, region, regulation, sep='_'),'.fa',sep=''),paste(paste('Regulated',geneVecName, region, sep='_'),'.fa',sep='')))
  seqinr::write.fasta(sequences=as.list(as.character(controlSeq)),names=names(controlSeq),file.out=ifelse(is.null(geneVec),paste(paste('Control_contrast',contrast, region, regulation, sep='_'),'.fa',sep=''),paste(paste('Control',geneVecName, region, sep='_'),'.fa',sep='')))
  #Run dreme
  #
  if(is.null(geneVec)){
    outdirTmp <- ifelse(is.null(stremeName),paste('stremeOut_contrast', contrast, region, regulation, sep='_'),stremeName)
  } else {
    outdirTmp <- ifelse(is.null(stremeName),paste('stremeOut', geneVecName, region, sep='_'),stremeName)
  }
  #
  streme_out <- memes::runStreme(input=ifelse(is.null(geneVec),paste(paste('Regulated_contrast',contrast, region, regulation, sep='_'),'.fa',sep=''),paste(paste('Regulated',geneVecName, region, sep='_'),'.fa',sep='')), control=ifelse(is.null(geneVec),paste(paste('Control_contrast',contrast, region, regulation, sep='_'),'.fa',sep=''),paste(paste('Control',geneVecName, region, sep='_'),'.fa',sep='')), alph="dna", outdir = outdirTmp,minw=minwidth)
  
  #Extract only these below threshold
  streme_out <- streme_out[streme_out$pval < stremeThreshold,]
  if(!is.null(tomtom_database)){
    if(nrow(streme_out)>0){
      if(is.null(geneVec)){
        outdirTmp2 <- ifelse(is.null(tomtomName),paste('tomtomOut_contrast', contrast, region, regulation, sep='_'),tomtomName)
      } else {
        outdirTmp2 <- ifelse(is.null(tomtomName),paste('tomtomOut', geneVecName, region, sep='_'),tomtomName)
      }
      #Run tomtom
      tomtom_out <- memes::runTomTom(streme_out,database=tomtom_database, outdir = outdirTmp2)
    } else {
      tomtom_out <- streme_out
    }
    return(tomtom_out)
  } else {
    #
    return(streme_out)
  }
}
