anota2seqGetDirectedRegulations <- function(ads) {
  n_contrasts <- ncol(ads@contrasts)
  regModeList <- vector("list", length = n_contrasts)
  
  for (i in 1:n_contrasts) {
    rvm_data <- ads@selectedTranslation@selectedRvmData[[i]]
    buffering_data <- ads@selectedBuffering@selectedRvmData[[i]]
    abundance_data <- ads@mRNAAbundance@translatedmRNA[[i]]
    total_mrna_data <- ads@selectedTotalmRNA@selectedRvmData[[i]]
    
    translationUp <- rvm_data[rvm_data$apvEff > 0 & rvm_data$singleRegMode == "translation", ]
    translationDown <- rvm_data[rvm_data$apvEff < 0 & rvm_data$singleRegMode == "translation", ]
    translatedmRNAUp <- rvm_data[rvm_data$apvEff > 0, ]
    translatedmRNADown <- rvm_data[rvm_data$apvEff < 0, ]
    bufferingmRNAUp <- buffering_data[buffering_data$apvEff > 0 & buffering_data$singleRegMode == "buffering", ]
    bufferingmRNADown <- buffering_data[buffering_data$apvEff < 0 & buffering_data$singleRegMode == "buffering", ]
    mRNAAbundanceUp <- abundance_data[abundance_data$apvEff > 0 & abundance_data$singleRegMode == "abundance", ]
    mRNAAbundanceDown <- abundance_data[abundance_data$apvEff < 0 & abundance_data$singleRegMode == "abundance", ]
    totalmRNAUp <- total_mrna_data[total_mrna_data$apvEff > 0, ]
    totalmRNADown <- total_mrna_data[total_mrna_data$apvEff < 0, ]
    
    regModeList[[i]] <- list(
      "translationUp" = rownames(translationUp),
      "translationDown" = rownames(translationDown),
      "translatedmRNAUp" = rownames(translatedmRNAUp),
      "translatedmRNADown" = rownames(translatedmRNADown),
      "bufferingmRNAUp" = rownames(bufferingmRNAUp),
      "bufferingmRNADown" = rownames(bufferingmRNADown),
      "mRNAAbundanceUp" = rownames(mRNAAbundanceUp),
      "mRNAAbundanceDown" = rownames(mRNAAbundanceDown),
      "totalmRNAUp" = rownames(totalmRNAUp),
      "totalmRNADown" = rownames(totalmRNADown)
    )
  }
  return(regModeList)
}


###
checkAvailableVersions <- function(species) {
  # Define the base directory for annotation files
  base_dir <- system.file("extdata/annotation/refseq", package = "anota2seqUtils")
  
  # List existing species
  curr_tmp <- list.files(base_dir)
  
  # Check if the input species is a valid species name
  if (!species %in% curr_tmp) {
    stop("Invalid species name. The pre-prepared annotation file for that species does not exist. Please use a valid species name.")
  } else {
    # Create the directory path for the specified species
    species_dir <- file.path(base_dir, species)
    
    # List available versions for the specified species
    list_versions <- list.files(path = species_dir)
    
    if (length(list_versions) == 0) {
      cat("No annotation versions found for", species, "\n")
    } else {
      cat("Available versions for", species, ":\n")
      print(list_versions)
    }
  }
}

###For gff 
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

gffRead <- function(gffFile, nrows = -1){
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  #stopifnot(!anyis.na(gff$start)), !anyis.na(gff$end)))
  return(gff)
}

extGff <- function(gff){
  #
  gff <- gff[gff$feature=='mRNA' | gff$feature=='transcript',]
  gff$transID_ver <- getAttributeField(gff$attributes, "transcript_id")
  gff$geneID <- getAttributeField(gff$attributes, "gene")
  gff <- with(gff,cbind(gff,reshape2::colsplit(gff$transID_ver,pattern="\\.",names = c('transID','version'))))
  bed <- data.frame(id=gff$transID,chr=gff$seqname, strand=gff$strand, start=gff$start, end=gff$end, transVer=gff$version, geneID=gff$geneID)
  #Extract only protein coding
  bed <- bed[grepl('NM_',bed$id),]
  bed <- bed[grepl('NC_',bed$chr),]
  bed <- subset(bed, !duplicated(id))
  #
  return(bed)
}


gSel <- function(annot,ads,customBg,geneList){
  if(!is.null(ads)){
    bg <- row.names(ads@dataP)
    #if(!is.null(geneList)){
    #  bg <- unique(c(bg, as.character(unlist(geneList))))
    #}
    annotOut <- annot[annot$geneID %in% bg,]
  } else {
    if(!is.null(customBg)){
      #add here to be sure that all genelist are in bg
      bg <- customBg
      if(!is.null(geneList)){
        bg <- unique(c(bg, as.character(unlist(geneList))))
      }
      annotOut <- annot[annot$geneID %in% bg,]
    } else {
      annotOut <- annot[annot$geneID %in% as.character(unlist(geneList)),]
    }
  }
  return(annotOut)
}

regSel <- function(annot, region, ext=FALSE){
  nc <- grep(region, colnames(annot))
  #
  seqTmp <- annot[,nc]
  lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
  #
  if(isTRUE(ext)){
    seq <- list()
    seq[[1]] <- annot$CDS_seq
    seq[[2]] <- annot$UTR3_seq
    extSeq <- combSeq(seqIn = seq)
    extSeq <- unlist(extSeq)
    #
    annotOut <- cbind(annot[,c(1:2)], seqTmp,lenTmp,extSeq)
  } else {
    annotOut <- cbind(annot[,c(1:2)], seqTmp,lenTmp)
  }
  #
  return(annotOut)
}

isoSel <- function(annot, method){
  #Select per gene level
  if(method=='shortest'){
    annotOut <- as.data.frame(annot %>% group_by(geneID) %>% dplyr::slice(which.min(lenTmp)))
  } else if(method=='longest'){
    annotOut <- as.data.frame(annot %>% group_by(geneID) %>% dplyr::slice(which.max(lenTmp)))
  } else {
    annotOut <- as.data.frame(annot %>% group_by(geneID) %>% dplyr::slice_sample(n = 1))
  }
  return(annotOut)
}

resSel <- function(vIn, ads=NULL, regulation=NULL, regulationGen=NULL, contrast, customBg, geneList){
  resOut <- list()
  #Extract all results
  if(!is.null(ads)){
    results <- anota2seqGetDirectedRegulations(ads)
    #
    if(!is.null(regulationGen)){
      regulation <- names(results[[contrast]])[grepl(regulationGen, names(results[[contrast]]))]
      contrast <- rep(contrast,length(regulation))
      res <- vector("list", length = length(regulation))
      for(i in unique(contrast)){
        resTmp <- results[[i]][regulation[contrast==i]]
        res[which(contrast==i)] <- resTmp
      }
      names(res) <- paste(regulation, paste('c', contrast,sep=''), sep='_')
    } 
    
    if (!is.null(regulation)){
      res <- vector("list", length = length(regulation))
      for(i in unique(contrast)){
        resTmp <- results[[i]][regulation[contrast==i]]
        res[which(contrast==i)] <- resTmp
      }
      names(res) <- paste(regulation, paste('c', contrast,sep=''), sep='_')
    } else {
      res <- list()
      for(i in 1:length(results)){
        resTmp <- results[[i]]
        names(resTmp) <- paste(names(resTmp), paste('c', i,sep=''), sep="_")
        res <- append(res, resTmp)
      }
    } 
    #if(!is.null(geneList)){
    #  res <- append(res, geneList)
    #}
    #if(isTRUE(addBg)){
    resOut[[1]] <- vIn
    #}
    for(i in 1:length(res)){
      resOut[[names(res)[i]]] <- vIn[names(vIn) %in% res[[i]]]
    }
    #if(isTRUE(addBg)){
    names(resOut)[1] <- 'background'
    #}
  } else {
    res <- geneList
    if(!is.null(customBg)){
      resOut[[1]] <- vIn
      for(i in 1:length(res)){
        resOut[[names(res)[i]]] <- vIn[names(vIn) %in% res[[i]]]
      }
      names(resOut)[1] <- 'background'
    } else {
      for(i in 1:length(res)){
        resOut[[names(res)[i]]] <- vIn[names(vIn) %in% res[[i]]]
      }
    }
  }
  resOut <- resOut[lapply(resOut,length)>0]
  return(resOut)
}

coloursSel <- function(resOut, geneList, geneListcolours){
  coloursOut <- as.character()
  if(is.null(geneList)){
    AnotaColours <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],RColorBrewer::brewer.pal(8,"Reds")[c(2,6)],RColorBrewer::brewer.pal(8,"Greens")[c(4,8)], RColorBrewer::brewer.pal(8,"Greens")[c(2,6)],RColorBrewer::brewer.pal(8,"Blues")[c(4,8)])
    names(AnotaColours) <- c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown","bufferingmRNAUp","bufferingmRNADown")
    #
    coloursOut <- AnotaColours[gsub("\\_.*","",names(resOut))]
    coloursOut[1] <- 'grey75'
  } else {
    if(names(resOut)[1] == 'background'){
      coloursOut <- c('grey65', geneListcolours)
    } else {
      coloursOut <- geneListcolours
    }
  } 
  return(coloursOut)
}

addStats <- function(comparisons, plotType, resOut, coloursOut){
  #
  if(plotType == "ecdf"){
    tableOut <- matrix(NA, nrow = length(comparisons), ncol = 5)
    colnames(tableOut) <- c("signature", "Wilcox_pval", "q25", "q50", "q75")
    colOut <- as.character()
  }
  for (j in 1:length(comparisons)) {
    if (names(resOut)[1] == 'background') {
      compTmp <- comparisons[[j]] + 1
    } else {
      compTmp <- comparisons[[j]]
    }
    # stats
    pvalTmp <- format(as.numeric(wilcox.test(resOut[[compTmp[1]]], resOut[[compTmp[2]]], alternative = "two.sided")[3]), scientific = TRUE, digits = 2)
    #
    if (plotType == "boxplot" | plotType == "violin") {
      yposTmp <- ifelse(range(as.numeric(unlist(resOut)))[2] <= 1, 1.1,range(as.numeric(unlist(resOut)))[2] + j*5)
      rect(xleft = compTmp[1], xright = compTmp[2], ybottom = yposTmp, ytop = yposTmp, lwd = 2)
      #
      text(sum(compTmp) / 2, ifelse(range(as.numeric(unlist(resOut)))[2] <= 1,yposTmp + 0.05,yposTmp + 2.5), pvalTmp, cex = 0.75)
    } else if (plotType == "ecdf") {
      tableOut[j, 1] <- paste(names(resOut)[compTmp[2]], "vs", names(resOut)[compTmp[1]], sep = " ")
      tableOut[j, 2] <- pvalTmp
      
      # Calculate percentiles
      tmpBg <- sort(resOut[[compTmp[1]]])
      ecdfBg <- 1:length(tmpBg) / length(tmpBg)
      bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
      bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
      bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
      
      # Calculate percentiles for second and difference from background
      tmpSign <- sort(resOut[[compTmp[2]]])
      ecdfSign <- 1:length(tmpSign) / length(tmpSign)
      tableOut[j, 3] <- format(tmpSign[which(ecdfSign >= 0.25)[1]] - bg_025, digits = 2)
      tableOut[j, 4] <- format(tmpSign[which(ecdfSign >= 0.5)[1]] - bg_05, digits = 2)
      tableOut[j, 5] <- format(tmpSign[which(ecdfSign >= 0.75)[1]] - bg_075, digits = 2)
      #
      if (length(which(grepl("background", c(names(resOut)[compTmp[2]], names(resOut)[compTmp[1]])))) > 0) {
        colT <- gsub("\\_.*", "", names(resOut)[compTmp][which(names(resOut)[compTmp] != "background")])
        colOut[j] <- coloursOut[colT]
      } else {
        colOut[j]  <- "white"
      }
      xlim_min <- floor(quantile(as.numeric(unlist(resOut)), 0.01))
    }
  }
  if (plotType == "ecdf"){
    plotrix::addtable2plot(xlim_min, 1.01, tableOut, bty = "n", display.rownames = FALSE, hlines = FALSE, vlines = TRUE, title = "", cex = 0.7, bg = colOut, xpad = 0.1, ypad = 1.4, xjust = 0, yjust = 1)
  }
}


# 
plotBoxplots <- function(resOut, coloursOut, comparisons) {
  par(mar = c(8, 12, 12, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
  xlimIn <- c(0.5, length(resOut) + 1.5)
  
  plot(1, 1, xlim = xlimIn, ylim = c(0, range(as.numeric(unlist(resOut)))[2] + (1.25 * length(comparisons))), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
  
  axis(side = 2, font = 2, las = 2, lwd = 2, at = sapply(c(1, 25, 100, 200, 400, 1000, 4000, 25000), log2), labels = c(0, 25, 100, 200, 400, 1000, 4000, 25000))
  mtext(side = 2, line = 6, "Log2 length", col = "black", font = 2, cex = 1.7, at = median(as.numeric(unlist(resOut))))
  text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOut), xpd = NA, cex = 0.9, srt = 45, adj = 1)
  
  if (names(resOut)[1] == 'background') {
    abline(lty = 5, h = median(resOut[[1]]))
  }
  #
  for (i in 1:length(resOut)) {
    boxplot(resOut[[i]], add = TRUE, at = i, col = coloursOut[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
    text(i, 0, round(mean(antilog(resOut[[i]], 2), 0)), font = 2)
  }
  if(!is.null(comparisons)){
    addStats(comparisons, plotType='boxplot', resOut, coloursOut)
  }
}

plotViolin <- function(resOut, coloursOut, comparisons) {
  par(mar = c(8, 12, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
  xlimIn <- c(0.5, length(resOut) + 1.5)
  
  plot(1, 1, xlim = xlimIn, ylim = c(0, range(as.numeric(unlist(resOut)))[2] + (1.25 * length(comparisons))), xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
  
  axis(side = 2, font = 2, las = 2, lwd = 2, at = sapply(c(1, 25, 100, 200, 400, 1000, 4000, 25000), log2), labels = c(0, 25, 100, 200, 400, 1000, 4000, 25000))
  mtext(side = 2, line = 6, "Log2 length", col = "black", font = 2, cex = 1.7, at = median(as.numeric(unlist(resOut))))
  text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOut), xpd = NA, cex = 0.9, srt = 45, adj = 1)
  
  if (names(resOut)[1] == 'background') {
    abline(lty = 5, h = median(resOut[[1]]))
  }
  #
  for (i in 1:length(resOut)) {
    vioplot::vioplot(resOut[[i]], add = TRUE, at = i, col = coloursOut[i], xaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
    text(i, 0, round(mean(antilog(resOut[[i]], 2), 0)), font = 2)
  }
  if(!is.null(comparisons)){
    addStats(comparisons, plotType='violin', resOut, coloursOut)
  }
}

# Helper function for plotting ECDF
plotEcdf <- function(resOut, coloursOut, comparisons) {
  xlim_min <- as.numeric(quantile(as.numeric(unlist(resOut)), 0.01))
  xlim_max <- as.numeric(quantile(as.numeric(unlist(resOut)), 0.99))
  par(mar = c(5, 5, 8, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
  
  #
  plot(ecdf(resOut[[1]]), col = coloursOut[1], main = "", xlab = "", ylab = "", verticals = TRUE, do.p = FALSE, lwd = 3, bty = "n", yaxt = "none", font = 2, xlim = c(xlim_min, xlim_max), xaxt = "none")
  
  mtext(side = 1, line = 4, "Log2 length", col = "black", font = 2, cex = 1.2)
  mtext(side = 2, line = 3, "Fn(x)", col = "black", font = 2, cex = 1.2)
  
  axis(side = 1, seq(floor(xlim_min), ceiling(xlim_max), 1), font = 2, lwd = 2)
  axis(side = 2, seq(0, 1, 0.2), font = 2, las = 2, lwd = 2)
  for (i in 2:length(resOut)) {
    lines(ecdf(resOut[[i]]), col = coloursOut[i], main = "", xlab = "", verticals = TRUE, do.p = FALSE, lwd = 4)
  }
  if(!is.null(comparisons)){
    addStats(comparisons, plotType='ecdf', resOut, coloursOut)
  }
}



antilog<-function(lx,base){ 
  lbx<-lx/log(exp(1),base=base) 
  result<-exp(lbx) 
  result 
} 


#calculate desired content
calc_content <- function(x, nuc){
  contTmp <- stringr::str_count(x, nuc)
  contentOut <- contTmp / stringr::str_length(x) * 100
  return(contentOut)
}

####
subset_seq <- function(x, pos, subregionSel){
  #
  if(length(seqinr::s2c(x)) >= abs(pos)){
    #
    if(pos>0){
      if(subregionSel=='select'){
        seqTmp <- seqinr::c2s(seqinr::s2c(x)[1:pos])
      } else if(subregionSel=='exclude'){
        seqTmp <- seqinr::c2s(seqinr::s2c(x)[(pos+1):length(seqinr::s2c(x))])
      }
    } else if (pos<0){
      if(subregionSel=='select'){
        seqTmp <- seqinr::c2s(seqinr::s2c(x)[(length(seqinr::s2c(x))-abs(pos)+1):length(seqinr::s2c(x))])
      } else if(subregionSel=='exclude'){
        seqTmp <- seqinr::c2s(seqinr::s2c(x)[1:(length(seqinr::s2c(x))-abs(pos))])
      }
    }
  } else {
    seqTmp <- NA
  }
  return(seqTmp)
}

#
calc_motif <- function(x, motifIn, dist, unit){
  seqTmp <- x
  #
  lenTmp <- motifLenCalc(motifIn)
  len <- ifelse(dist == 1, lenTmp, lenTmp + (dist - 1))
  #
  motOut  <- seqinr::words.pos(motifIn,seqTmp)
  if(length(motOut)>0){
    #Check overlapping and collapse them
    gROut <- GenomicRanges::reduce(GenomicRanges::GRanges(seqnames='tmp', ranges=IRanges::IRanges(start=motOut,end=motOut)),min.gapwidth=len)
    #
    if(unit == 'number'){
      nMot <- length(gROut)
    } else if (unit == 'position') {
      nMot <- list()
      nMot[["start"]] <- as.numeric(start(gROut@ranges))
      nMot[["end"]] <- as.numeric(end(gROut@ranges)) + lenTmp - 1
    }
  } else {
    nMot <- ifelse(unit == "number", 0, NA)
  }
  return(nMot)
}

#Combine 
calc_g4 <- function(x,min_score){
  seqTmp <- DNAString(x)
  predTmp <- pqsfinder::pqsfinder(seqTmp, min_score = min_score, strand = '+', verbose=F)
  nMot <- nrow(predTmp@elementMetadata)
  return(nMot)
}


#convert IUPAC code
convertIUPAC <- function(motif){
  #
  tmpConv <- toupper(motif)
  #Improve one day
  tmpConv <- gsub('N','[ACGT]',tmpConv)
  tmpConv <- gsub('U','T',tmpConv)
  tmpConv <- gsub('R','[AG]',tmpConv)
  tmpConv <- gsub('Y','[CT]',tmpConv)
  tmpConv <- gsub('K','[GT]',tmpConv)
  tmpConv <- gsub('M','[AC]',tmpConv)
  tmpConv <- gsub('S','[CG]',tmpConv)
  tmpConv <- gsub('W','[AT]',tmpConv)
  tmpConv <- gsub('B','[CGT]',tmpConv)
  tmpConv <- gsub('D','[GAT]',tmpConv)
  tmpConv <- gsub('H','[ACT]',tmpConv)
  tmpConv <- gsub('V','[ACG]',tmpConv)
  #
  return(tmpConv)
}

motifLenCalc <- function(motif){
  lenTmp <- length(seqinr::s2c(gsub("\\[|\\]", "", motif)))
  #
  bracSel <- unlist(regmatches(motif, gregexpr("(?<=\\[).*?(?=\\])", motif, perl=T)))
  lenB <- length(seqinr::s2c(seqinr::c2s(bracSel)))
  #
  nNotInBrackets <- lenTmp  - lenB
  #
  lenOut <- length(bracSel) + nNotInBrackets
  return(lenOut)
}


#prepprotein
replaceProtAmbig <- function(motif){
  #
  tmpConv <- toupper(motif)
  #Improve one day
  tmpConv <- gsub('X','[ACDEFGHIKLMNPQRSTVWY]',tmpConv)
  #
  return(tmpConv)
}

#
codonCount <- function(seq, gene, codonN=1){
  #
  if(codonN==1){
    #
    tmpEff <- seqinr::uco(seqinr::s2c(seq),index = "eff")
    tmpFreq <- seqinr::uco(seqinr::s2c(seq),index = "freq")
    #
    tmpCodon <- data.frame(geneID=gene,codon=toupper(names(tmpEff)),AA=seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(toupper(names(tmpEff)))))), codonCount=as.numeric(tmpEff),codonFreq=as.numeric(tmpFreq))
    tmpCodon <- tmpCodon %>% group_by(AA) %>% mutate(AACountPerGene=sum(codonCount))
  } else if(codonN > 1){
    seqIn <- seqinr::s2c(tolower(seq))
    
    tmpEff <- seqinr::count(seq = seqIn, wordsize = 3*codonN, start = 3, by = 3,freq = FALSE)
    tmpFreq<- length(seqIn)/(3*codonN)
    #tmpFreq <- seqinr::count(seq = seqIn, wordsize = 3*codonN, start = 3, by = 3,freq = TRUE)
    #
    #indSel <- which(tmpEff>0)
    #
    #tmpCodon <- data.frame(geneID=gene,codon=toupper(names(tmpEff)[indSel]),AA=NA, codonCount=as.numeric(tmpEff)[indSel],codonFreq=as.numeric(tmpFreq)[indSel])
    tmpCodon <- data.frame(geneID=gene,codon=toupper(names(tmpEff)),AA=NA, codonCount=as.numeric(tmpEff),codonFreq=as.numeric(tmpFreq))
    
  }
  #
  return(tmpCodon)
}

#
statOnDf <- function(df, # dataframe with summed codon counts for each regulation
                     regs, # regulations you want to be tested (only 2 i.e. trans Up vs trans Down)
                     analysis #codon or AA
){
  #
  oddRatioOut <- list()
  #
  uniqAA <- unique(as.character(df$AA))
  
  if(analysis=='codon'){
    #df$p.value <- 1
    ##
    ## Create a list that holds all the 2*2 tables on which the fisher test is performed
    fisherList <- rep(list(NA),length(unique(df$codon)))
    testList <- rep(list(NA),length(unique(df$codon)))
    names(testList) <- names(fisherList)<- unique(df$codon)
  
    for(AAind in 1:length(uniqAA)){
      #
      tmpDf <- df[df$AA == uniqAA[AAind],]
      #
      codons <- as.character(df$codon[df$AA == uniqAA[AAind]])
    
      # Check if the AA has 2 codons
      if(length(codons) > 1){
        for(cod in 1:length(codons)){
          codTmp <- codons[cod]
          c_codTmp <- tmpDf[tmpDf$codon==codTmp,regs]
          #calculate all that are not in that one
          r_codTmp <- apply(tmpDf[!tmpDf$codon==codTmp,regs],2,sum)
          #
          fisherIn <- rbind(c_codTmp,r_codTmp)
        
          ## do the fisher test on 2 codons
          fisherOut <- fisher.test(fisherIn)
          #
          oddRatioOut[[codTmp]] <- as.numeric(fisherOut$estimate)
        }
      } else {
        oddRatioOut[[codTmp]] <- NA
      }
    }
  } else if(analysis=='AA'){
    fisherList <- rep(list(NA),length(unique(df$AA)))
    testList <- rep(list(NA),length(unique(df$AA)))
    names(testList) <- names(fisherList)<- unique(df$AA)
    
    for(AAi in 1:length(uniqAA)){
      AATmp <- uniqAA[AAi]
      c_codTmp <- df[df$AA == AATmp,regs]
      #calculate all that are not in that one
      r_codTmp <- colSums(df[df$AA != AATmp,regs])
      
      fisherIn <- rbind(c_codTmp,r_codTmp)
      ## do the fisher test on 2 codons
      fisherOut <- fisher.test(fisherIn)
      
      #
      oddRatioOut[[AATmp]] <- as.numeric(fisherOut$estimate)
    }
  }
  oddRatioOut<- unlist(oddRatioOut)
  #
  return(oddRatioOut)
}


roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
     if(length(x) != 1) stop("'x' must be of length 1")
     10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

combSeq <- function(seqIn){
  seqTmp <- lapply(seqIn, function(x) lapply(x, function(y) seqinr::s2c(y)))
  #
  seqC <- do.call(Map, c(c, seqTmp))
  #
  seqOut <- lapply(seqC, function(x) seqinr::c2s(x))
  #
  return(seqOut)
}

calc_uORF <- function(seqTmp, ext, context , unit){
  #
  nTmpStart <- as.numeric()
  nTmpStop <- as.numeric()
  #
  #Look for ATG in context in 5UTR
  startOut  <- seqinr::words.pos(context, seqTmp) + 3
  if(length(startOut)>0){
    #Look for stop in 5UTR or full transcript
    if(!is.null(ext)){
      seqTmp2 <- seqinr::c2s(c(seqinr::s2c(seqTmp), seqinr::s2c(ext)))
    } else {
      seqTmp2 <- seqTmp
    }
    #
    stopOut <- seqinr::words.pos('TAA|TAG|TGA', seqTmp2)
    #
    #stopOut <- stopOut + 2
    if(length(stopOut)>0){
      #
      for(i in 1:length(startOut)){
        #Remove if stop is before start 
        stopTmp <- stopOut[which((stopOut - startOut[i])>0)]
        #
        potORF <- stopTmp - startOut[i]
        #check in frame 
        inFrameCheck <- potORF %% 3
        #Take first stop in frame
        stopTmp <- stopTmp[which(inFrameCheck==0)]
        #if maybe postition needed
        #if(length(stopOut)>0){
        #  stopOut <- min(stopOut)+2
        #}
        #
        if(length(stopTmp)>0){
          nTmpStart[i] <- startOut[i]
          nTmpStop[i] <- min(stopTmp)
        }
      }
      #
      if(unit == 'number'){
        nOut <- length(nTmpStart)
      } else if (unit == 'position') {
        nOut <- list() 
        nOut[["start"]] <- nTmpStart
        nOut[["end"]] <- nTmpStop
      }
    } else {
      nOut <- ifelse(unit == "number", 0, NA)
    }
  } else {
    nOut <- ifelse(unit == "number", 0, NA)
  }
  return(nOut)
}


convertSymbolToEntrezID <- function(geneList,
                                    species
){
  if(species=="human"){
    symbol2id <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL2EG)
  } else if(species=="mouse"){
    symbol2id <- as.list(org.Mm.eg.db::org.Mm.egSYMBOL2EG)
  }
  #          
  identifiersWithNoEntrezID <- c()
  
  symbols <- geneList
  if(length(symbols)>0){
    geneListEntrezID <- unique(na.omit(as.character(unlist(symbol2id[intersect(names(symbol2id), unique(symbols))]))))
    #identifiersWithNoEntrezID <- c(identifiersWithNoEntrezID, setdiff(unique(symbols), names(symbol2id)))
  }
  
  # if some element of the list have become empty (empty gene signature) because none of the symbols have an entrezid
  
  if(length(geneListEntrezID) == 0){
    stop("None of the provided gene symbols could be converted to entrezid. Please check the gene identifiers")
  }
  
  return(geneListEntrezID)
}

convertEntrezIDToSymbol <- function(entrezIDList,
                                    species
){
  if(species=="human"){
    id2symb <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
  } else if(species=="mouse"){
    id2symb <- as.list(org.Mm.eg.db::org.Mm.egSYMBOL)
  }
  #          
  identifiersWithNoEntrezID <- c()
  
  entrezIDs <- entrezIDList
  if(length(entrezIDs)>0){
    geneList <- unique(na.omit(as.character(unlist(id2symb[intersect(names(id2symb), unique(entrezIDs))]))))
    #identifiersWithNoEntrezID <- c(identifiersWithNoEntrezID, setdiff(unique(symbols), names(symbol2id)))
  }
  
  # if some element of the list have become empty (empty gene signature) because none of the symbols have an entrezid
  
  if(length(geneList) == 0){
    stop("None of the provided entrezIDs could be converted to symbols. Please check the identifiers")
  }
  
  return(geneList)
}


writeExcel <- function(listOfData=NULL, listNames=NULL, fileName=NULL){
  
  listToWrite <- list()
  for(i in 1:length(listOfData)){
    if(!is.null(listOfData[[i]])){
      listToWrite[[i]]<- as.data.frame(listOfData[[i]])
    }
    else{
      listToWrite[[i]] <- as.data.frame("none")
    }
  }
  names(listToWrite)<- listNames
  WriteXLS::WriteXLS(listToWrite,fileName,SheetNames =listNames, row.names=TRUE)
}

#For nodes sizes:
rescale <- function(x,a,b,c,d){
  c + (x-a)/(b-a)*(d-c)
}

#If the names are long, wrap them up
wrapNames <- function(s,w){
  as.character(sapply(s, FUN=function(x){
    paste(strwrap(x, width=w), collapse="\n")
  }))
}


layoutCalc <- function(Gobject, n){  
  Gtmp <- Gobject
  igraph::E(Gtmp)$weight <- 1
  
  attr <- cbind(id=1:igraph::vcount(Gtmp), val=n)
  Gtmp <- Gtmp + igraph::vertices(unique(attr[,2])) + igraph::edges(unlist(t(attr)), weight=0.25)
  
  lOut <- igraph::layout_nicely(Gtmp, weights=igraph::E(Gtmp)$weight)[1:igraph::vcount(Gobject),]
  return(lOut)
}

extractRegSeq <- function(annotSeq){
  #
  annotOut <- annotSeq
  #
  #5UTR
  seqSel <- as.character()
  for(i in 1:nrow(annotSeq)){
    seqTmp <- annotSeq$seq[i]
    #
    utr5len <- as.numeric(annotSeq$UTR5_len[i])
    #
    seqSel[i] <- seqinr::c2s(seqinr::s2c(seqTmp)[1:(utr5len)])
  }
  annotOut$UTR5_seq <- seqSel
  
  #CDS
  seqSel <- as.character()
  for(i in 1:nrow(annotSeq)){
    seqTmp <- annotSeq$seq[i]
    #
    cdsStart <- as.numeric(annotSeq$UTR5_len[i])+1
    cdsStop <- as.numeric(annotSeq$CDS_stop)[i]
    #
    seqSel[i] <- seqinr::c2s(seqinr::s2c(seqTmp)[cdsStart:cdsStop])
  }
  annotOut$CDS_seq <- seqSel
  
  #3UTR
  seqSel <- as.character()
  for(i in 1:nrow(annotSeq)){
    seqTmp <- annotSeq$seq[i]
    #
    cdsStop <- as.numeric(annotSeq$CDS_stop)[i]
    totalLength <- as.numeric(annotSeq$Total_len)[i]
    #
    seqSel[i] <- seqinr::c2s(seqinr::s2c(seqTmp)[(cdsStop+1):totalLength])
  }
  annotOut$UTR3_seq <- seqSel
  #
  return(annotOut)
}


#Convert to dpn
generateOut <- function(x, tmpList){
  if(length(tmpList[[x]])>0){
    tmpPos <- paste(names(tmpList[[x]]), collapse="\t")
    tmpData <- paste(as.numeric(tmpList[[x]]), collapse="\t")
  }
  tmpString <- paste(c(x, "ZZZ", tmpData, "ZZZ", tmpPos), collapse="\t")
  return(tmpString)
}

codPlot <- function(rust, name){
  rustOut <- matrix(NA, nrow=nrow(rust)-1, ncol=60)
  for(i in 1:nrow(rust[-62,])){
    obs <- as.numeric(rust[i,3:62])
    exp <- as.numeric(rust[i,2])
    ratio <- log2(obs/exp)
    rustOut[i,] <- ratio 
  }
  pdf(paste(name,'codons.pdf',sep='_'),height=6,width=6)
  par(mar=c(5,5,5,4),bty='l',font=2, font.axis=2, font.lab=2, cex.axis=0.8,cex.main=0.8,cex.lab=1)
  plot(seq(-40,19,1),rustOut[1,],type='l',xlim=c(-40,20),ylim=c(-3,3),xaxt='n',xlab='',ylab='Codon RUST ratio',col='black',lwd=2)
  axis(seq(-40,20,5),side=1,at=seq(-40,20,5),tick=F)
  for(i in 2:nrow(rustOut)){
    lines(seq(-40,19,1),as.numeric(rustOut[i,]),col='black',lwd=2)
  }
  lines(seq(-40,19,1),as.numeric(rust[62,3:62]),col='dodgerblue',lwd=2)
  abline(v=0,col='red')
  dev.off()
}


readRiboDpn <- function(dpnFile, dpn_path=NULL, annot, cds_filt=TRUE){
  cat("Reading and processing", dpnFile, "\n")
  
  ##initate output object which is a list for refseqs.
  dataList <- list()
  #
  tmpList <- scan(paste(dpn_path,dpnFile,sep='/'), what="", sep="\n", quiet=TRUE)
  tmpList <- strsplit(tmpList, "\tZZZ\t")
  names(tmpList) <- sapply(tmpList, function(x) x[1])
  refList <- lapply(tmpList, modList)
  #subset to per gene
  refList <- refList[names(refList) %in% annot$id]
  #
  if(isTRUE(cds_filt)){
    #
    tmpAllRefs <- names(refList)
    #
    cds_start <- as.numeric(sapply(annot$UTR5_seq, function(x) length(seqinr::s2c(x))))+1
    cds_end <- cds_start + annot$lenTmp #as.numeric(sapply(annot$CDS_seq, function(x) length(seqinr::s2c(x))))
    names(cds_start) <- names(cds_end) <- annot$id
    #
    refListMod <- sapply(tmpAllRefs, extractInCDS, tmpList=refList, cds_start=cds_start, cds_end=cds_end)
  } else {
    refListMod <- refList
  }
  ##some refseqs may then loose all the data. These are removed
  tmpSel <- unlist(lapply(refListMod, length))
  refListMod <- refListMod[tmpSel>0]
  
  ##Replace transIDs with geneIDs
  names(refListMod) <- annot$geneID[match(names(refListMod), annot$id)]
  
  cat("done\n")
  
  ##add summary to output and return
  dataList <- refListMod
  return(dataList)
}

#
modList <- function(x){
  data <- as.numeric(unlist(strsplit(x[2], "\t")))
  pos <-  unlist(strsplit(x[3], "\t"))
  names(data) <- as.numeric(pos)
  return(data)
}

extractInCDS <- function(x, tmpList, cds_start, cds_end){
  #
  stTmp <- as.numeric(cds_start[x])
  endTmp <- as.numeric(cds_start[x] + cds_end[x])
  #
  dataTmp <- tmpList[[x]]
  #
  dataTmp <- dataTmp[as.numeric(names(dataTmp)) > stTmp & as.numeric(names(dataTmp)) <= endTmp]
  #Calculate relative to start of CDS
  names(dataTmp) <- as.numeric(names(dataTmp)) - stTmp
  #
  return(dataTmp)
}

countRiboCodon <- function(dataList, annot){
  #
  codonList <- list()
  #
  for(i in 1:length(dataList)){
    dataTmp <- dataList[[i]]
    #
    codonsOut <- list()
    for(j in 1:length(dataTmp)){
      #
      id <- names(dataTmp[j])
      datT <- dataTmp[[j]]
      #extract only in frame
      datT <- datT[as.numeric(names(datT)) %% 3 == 0]
      #
      if(length(datT) > 0){
        seqT <- annotSel$CDS_seq[annotSel$geneID == id]
        #
        codonTmp <- as.numeric(datT)
        names(codonTmp) <- sapply(as.numeric(names(datT)), extract_seq, seqs=seqT)
        #
        tmpFreq <- as.numeric(codonTmp)/sum(as.numeric(codonTmp))
        
        tmpCodon <- data.frame(geneID=id,codon=toupper(names(codonTmp)),AA=seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(toupper(names(codonTmp)))))), codonCount=as.numeric(codonTmp),codonFreq=as.numeric(tmpFreq))
        tmpCodon <- tmpCodon %>% group_by(AA) %>% mutate(AACountPerGene=sum(codonCount))
        
        codonsOut[[id]] <- tmpCodon
      } 
      #else {
      #  codonTmp <- NA
      #  names(codonTmp) <- NA
      #}
      #codonsOut[[id]] <- codonTmp
    }
    codonList[[i]] <- codonsOut
  }
}
#      cdsS <- cdsLengths_start[transID]
#      cdsE <- cdsLengths_end[transID]
#      #
#      if(!is.null(seqT)){
#        #
#        transTmp <- dataTmp[[t]]
#        #Filter only these in 0 frame
#        #
#        transTmp <- transTmp[(as.numeric(names(transTmp))-as.numeric(cdsS)) %% 3 == 0]
#        codTmp <- sapply(names(transTmp), extract_seq)
#        names(codTmp) <- as.numeric(transTmp)
#      }
#      codObs <- as.numeric()
#      for(cod in 1:length(codFreq)){
#        codSum <- sum(as.numeric(names(codTmp[codTmp==names(codFreq[cod])])))
#        codObs[cod] <- codSum
#      }
#      names(codObs) <- names(codFreq) 
#      #rOut <- resid(lm(codObs ~ codFreq))
#      #outTmp[[transID]] <- rOut
#      outTmp[[transID]] <- codObs
#    }
#    dataCodon_0f[[d]] <- outTmp
#  }
#}

#
extract_seq <- function(pos, seqs){
  #
  start <- pos+1
  seq <- seqinr::s2c(seqs)
  #
  if(!is.null(seq)){
    seqv   <-  seqinr::c2s(seq[start:(start+2)])
    return(seqv)
  } else {
    return(NA)
  }
}

runMfold <- function(fastaFile){
  #
  nameTmp <- gsub('.fa','',fastaFile) 
  #
  seqsToFold <- seqinr::read.fasta(fastaFile)
  #
  logFile = paste(nameTmp,"logFile.txt",sep='_')
  cat("", file=logFile, append=FALSE, sep = "\n")
  
  resFile = paste(nameTmp,'foldEnergy.txt',sep='_')
  cat("id\tfold_energy\tlength", file=resFile, append=FALSE, sep = "\n")
  
  pb1 <- txtProgressBar(min=1, max=length(seqsToFold), style=3)
  for(seqs in 1:length(seqsToFold)){    
    if(seqinr::c2s(seqsToFold[[seqs]]) != "na"){
      seqinr::write.fasta(seqsToFold[[seqs]], as.string = FALSE, names = attributes(seqsToFold[[seqs]])$name, 
                          file.out = "currentSeq.fa", open = "w")
      
      tryCatch({
        system("mfold SEQ=currentSeq.fa > currentSeq.out")
        cat(paste(c(attributes(seqsToFold[[seqs]])$name, read.table("currentSeq.fa_1.ct", nrows = 1, 
                                                                    stringsAsFactors = FALSE, 
                                                                    fileEncoding="UTF-8", 
                                                                    header = F)[, c("V4","V1")]), ## the identifier is also in currentSeq.fa_1.ct but often with encoding/formatting issues
                  collapse = "\t"), file = resFile, append = TRUE, sep = "\n")
      }, error = function(e){
        cat(paste("error", attributes(seqsToFold[[seqs]])$name), file=logFile, append=TRUE, sep = "\n")
      })
      
    } else { ##seq is NA
      cat(paste(c(rep(NA, 2), attributes(seqsToFold[[seqs]])$name), collapse = "\t"), 
          file = resFile, append = TRUE, sep = "\n")
    }
    
    unlink(list.files(pattern = "currentSeq"))
    setTxtProgressBar(pb1, seqs)
  }
  close(pb1)
}

runLM <- function(dataIn, namesDf, allFeat, useCorel, nameOut, NetModelSel, isads=NULL, regulationGen=NULL){
  #
  #Univariate
  models_prerun <- lapply(colnames(dataIn)[-length(colnames(dataIn))], function(x) {
    anova(lm(substitute(effM ~ i, list(i = as.name(x))), data = dataIn))
  })
  #
  step1expl <- round(sapply(models_prerun, function(x) (x[1, 2] / (sum(x[1, 2], x[2, 2]))) * 100), 2)
  names(step1expl) <- namesDf$originalNames

  step1pval <- sapply(models_prerun, function(x) (x[1, 5]))
  names(step1pval) <- namesDf$originalNames
  
  step1pval_fdr <- p.adjust(step1pval)
  names(step1pval_fdr) <- namesDf$originalNames

  presel <- as.numeric(which(step1pval_fdr > 0.05 | is.na(step1pval_fdr)))

  if (length(presel) > 0) {
    step1expl <- step1expl[-presel]
    step1pval <- step1pval[-presel]
    step1pval_fdr <- step1pval_fdr[-presel]
  }

  if (!isTRUE(allFeat) & length(presel) > 0) {
    dataIn <- dataIn[, -presel]
    #featureNames <- featureNames[-presel]
  }
  #
  fval <- list()
  pval <- list()

  models <- lapply(colnames(dataIn)[-length(colnames(dataIn))], function(x) {
    anova(lm(substitute(effM ~ i, list(i = as.name(x))), data = dataIn))
  })
  #
  fvalTmp <- sapply(models, function(x) x[1, 4])
  names(fvalTmp) <- colnames(dataIn)[-length(colnames(dataIn))]
  fvalTmp[is.na(fvalTmp)] <- 0
  fval[[1]] <- fvalTmp

  pvalTmp <- sapply(models, function(x) x[1, 5])
  names(pvalTmp) <- colnames(dataIn)[-length(colnames(dataIn))]
  pvalTmp[is.na(pvalTmp)] <- 1
  pval[[1]] <- pvalTmp

  #
  bestTmp <- which.max(sapply(models, function(x) x[1, 4]))
  outTmp <- which(sapply(models, function(x) x[1, 5]) > 0.05 | is.na(sapply(models, function(x) x[1, 5])))

  bestSel <- colnames(dataIn)[-length(colnames(dataIn))][bestTmp]
  outSel <- colnames(dataIn)[-length(colnames(dataIn))][outTmp]

  i <- 1
  # 
  while (length(colnames(dataIn)[-length(colnames(dataIn))][!colnames(dataIn)[-length(colnames(dataIn))] %in% c(bestSel, outSel)]) > 0) {
    #
    i <- i + 1
    #
    tmpIn <- colnames(dataIn)[-length(colnames(dataIn))][!colnames(dataIn)[-length(colnames(dataIn))] %in% c(bestSel, outSel)]
    models2 <- list()
    #
    x_sel <- paste(bestSel, collapse = " + ")
  
    for (j in 1:(length(colnames(dataIn)[-length(colnames(dataIn))]) - length(c(bestSel, outSel)))) {
      varx <- colnames(dataIn)[-length(colnames(dataIn))][!colnames(dataIn)[-length(colnames(dataIn))] %in% c(bestSel, outSel)][j]
      design <- as.formula(paste(paste("effM", paste(x_sel, collapse = " + "), sep = "~"), varx, sep = "+"))
      models2[[j]] <- anova(lm(design, data = dataIn))
    }
    #
    #if(length(which(is.na(lm(design, data = dataIn)$coefficients)))>0){
    #  fvalTmp  <- NA
    #  fval[[i]] <- fvalTmp
    #} else {
      fvalTmp <- sapply(models2, function(x) x[nrow(x) - 1, 4])
      names(fvalTmp) <- colnames(dataIn)[-length(colnames(dataIn))][!colnames(dataIn)[-length(colnames(dataIn))] %in% c(bestSel, outSel)]
      fval[[i]] <- fvalTmp
    #}
  
    #
    #if(length(which(is.na(lm(design, data = dataIn)$coefficients)))>0){
    #  pvalTmp  <- NA
    #  pval[[i]] <- pvalTm <- fvalTmp
    #} else {
      pvalTmp <- sapply(models2, function(x) x[nrow(x) - 1, 5])
      names(pvalTmp) <- colnames(dataIn)[-length(colnames(dataIn))][!colnames(dataIn)[-length(colnames(dataIn))] %in% c(bestSel, outSel)]
      pval[[i]] <- pvalTmp
    #}
    
    #
    bestTmp <- which.max(sapply(models2, function(x) x[nrow(x) - 1, 4]))
    outTmp <- which(sapply(models2, function(x) x[nrow(x) - 1, 5]) > 0.05)
    #
    if (length(outTmp) > 0) {
      if (!bestTmp %in% outTmp) {
        bestSel <- append(bestSel, tmpIn[bestTmp])
      }
    } else {
      bestSel <- append(bestSel, tmpIn[bestTmp])
    }
    outSel <- append(outSel, tmpIn[outTmp])
    #outSel <- append(outSel, 'a8')
    
  }
  #
  tmp <- t(plyr::ldply(fval, rbind))
  tmp <- tmp[c(bestSel, rev(outSel)), ]

  tmp_pval <- t(plyr::ldply(pval, rbind))
  tmp_pval <- tmp_pval[c(bestSel, rev(outSel)), ]
  
  tb1Out <- round(apply(tmp, 2, as.numeric), 2)
  row.names(tb1Out) <- namesDf$originalNames[match(row.names(tmp), namesDf$newNames)]
  tb1Out[is.na(tb1Out)] <- ""
  colnames(tb1Out) <- paste("step", seq(1, ncol(tb1Out), 1), sep = "")

  colours <- matrix("white", nrow(tb1Out), ncol(tb1Out))
  linkIn <- matrix(NA, nrow(tb1Out), ncol(tb1Out))
  
  for (k in 2:nrow(colours)) {
    trow <- as.numeric(na.omit(as.numeric(tb1Out[k, ])))
    #
    tsel_tab <- which((lag(trow) / trow) > 2)
    tsel_net <- lag(trow) - trow
    #
    colours[k, tsel_tab] <- "#FDE0C5"
    linkIn[k, 1:length(tsel_net)] <- tsel_net
  }
  
  for (k in 1:ncol(colours)) {
    colours[k, k] <- "#B0F2BC"
  }
  #
  colours[which(tmp_pval > 0.05)] <- "#B14D8E"
  linkIn[which(abs(linkIn) < covarFilt)] <- NA
  row.names(linkIn) <- row.names(tmp)
  #
  tt1 <- gridExtra::ttheme_default(core = list(fg_params = list(fontface = c(rep("plain", ncol(tb1Out)))), bg_params = list(fill = colours, col = "black")))
  tg1 <- gridExtra::tableGrob(tb1Out, theme = tt1)
  
  #check whether any of the cooeficients is NA as it will be removed
  
  #
  
  varExpl <- anova(lm(as.formula(paste("effM", paste(bestSel, collapse = " + "), sep = "~")), data = dataIn))
  varExpldepend <- numeric()
  for (i in 1:length(bestSel)) {
    tmpFeatI <- bestSel[i]
    #
    varExpldepend[i] <- round((varExpl[i, 2] / sum(varExpl$`Sum Sq`)) * 100, 2)
  }
  names(varExpldepend) <- bestSel
  
  #
  varExplIndepend <- numeric()
  for (i in 1:length(bestSel)) {
    tmpFeat <- bestSel[i]
    restFeat <- bestSel[-i]
    
    design <- as.formula(paste(paste("effM", paste(restFeat, collapse = " + "), sep = "~"), tmpFeat, sep = "+"))
    tmpM <- anova(lm(design, data = dat))

    varExplIndepend[i] <- round((tmpM[nrow(tmpM) - 1, 2] / sum(tmpM$`Sum Sq`)) * 100, 2)
  }
  names(varExplIndepend) <- bestSel
  
  #
  varExplIndepend2 <- numeric()
  #
  if (length(presel) > 0) {
    step1sel <- namesDf$newNames[-presel]
  } else {
    step1sel <- namesDf$newNames
  }
  for (i in 1:length(step1sel)) {
    tmpFeat2 <- step1sel[i]
    restFeat2 <- step1sel[-i]
    
    design2 <- as.formula(paste(paste("effM", paste(restFeat2, collapse = " + "), sep = "~"), tmpFeat2, sep = "+"))
    tmpM2 <- anova(lm(design2, data = dat))

    varExplIndepend2[i] <- round((tmpM2[nrow(tmpM2) - 1, 2] / sum(tmpM2$`Sum Sq`)) * 100, 2)
  }
  names(varExplIndepend2) <- step1sel
  
  #
  tb2out <- data.frame(Features = namesDf$originalNames[match(bestSel, namesDf$newNames)], Pvalue = format(varExpl$`Pr(>F)`[1:length(bestSel)], scientific = T, digits = 2), VarianceExplained_Omnibus = as.numeric(varExpldepend), VarianceExplained_Adjusted = as.numeric(varExplIndepend))
  tg2 <- gridExtra::tableGrob(tb2out, rows = NULL)
  
  tb3out <- data.frame(Features = names(step1expl), Pvalue_Univariate = format(as.numeric(step1pval), scientific = T, digits = 2), FDRvalue_Univariate = format(as.numeric(step1pval_fdr), scientific = T, digits = 2), VarianceExplained_Univariate = as.numeric(step1expl))
  tb3out <- tb3out[with(tb3out, order(-tb3out$VarianceExplained_Univariate)), ]
  
  tg3 <- gridExtra::tableGrob(tb3out, rows = NULL)
  
  # 
  tb4out <- data.frame(Features = namesDf$originalNames[match(names(varExplIndepend2), namesDf$newNames)], VarianceExplained_IndependentAll = as.numeric(varExplIndepend2))
  tb4out <- tb4out[with(tb4out, order(-tb4out$VarianceExplained_IndependentAll)), ]
  tg4 <- gridExtra::tableGrob(tb4out, rows = NULL)
  
  pdf(paste(nameOut, "varexpl_independAll.pdf", sep = "_") , width = dim(tg4)[2] + dim(tg4)[2] + 4, height = dim(tg4)[1] / 2, useDingbats = F)
  gridExtra::grid.arrange(tg4, ncol = 1, nrow = 1)
  dev.off()
  
  # 
  linkOut <- list()
  for (i in 2:ncol(linkIn)) {
    tmpIn <- linkIn[, i]
    #
    tmpOut <- as.numeric(tmpIn)[which(!is.na(tmpIn))]
    if (length(tmpOut) > 0) {
      names(tmpOut) <- paste(row.names(linkIn)[i - 1], names(tmpIn)[which(tmpIn != 0)], sep = "_")
    }
    #
    linkOut[[i - 1]] <- tmpOut
  }
  linkOut <- as.data.frame(unlist(linkOut))
  linkOut <- with(linkOut, cbind(linkOut, reshape2::colsplit(row.names(linkOut), pattern = "\\_", names = c("from", "to"))))
  rownames(linkOut) <- NULL
  colnames(linkOut)[1] <- "weight"
  linkOut <- linkOut[, c(2, 3, 1)]
  
  #
  if(isTRUE(useCorel)){
    for(i in 1:nrow(linkOut)){
      f1tmp <- linkOut[i,]$from
      f2tmp <- linkOut[i,]$to
      #
      f1dat <- dat[,colnames(dat)==f1tmp]
      f2dat <- dat[,colnames(dat)==f2tmp]
      #
      corOut <- cor(f1dat,f2dat)
      linkOut$weight[i] <- corOut
    }
    ##
    tb5Out <- linkOut
    tb5Out$weight <- round(tb5Out$weight,2)
    tb5Out$from <- namesDf$originalNames[match(tb5Out$from, namesDf$newNames)]
    tb5Out$to <- namesDf$originalNames[match(tb5Out$to, namesDf$newNames)]
    colnames(tb5Out) <- c('featureFrom','featureTo','CorCoef')
    ##
    tg5 <- gridExtra::tableGrob(tb5Out, rows = NULL)
  
    pdf(paste(nameOut, "FinalModel.pdf", sep = "_"), width = dim(tg3)[2] + dim(tg1)[2] + dim(tg2)[2] + dim(tg5)[2] + 14.5, height = dim(tg1)[1] / 2, useDingbats = F)
    gridExtra::grid.arrange(tg3, tg1, tg2, tg5, ncol = 4, nrow = 1, padding = 0, top = 0, left = 0)
    grid::grid.text(paste("Total variance explained: ", sum(as.numeric(varExpldepend)), "%", sep = ""), x = grid::unit(0.75, "npc"), y = grid::unit(0.90, "npc"), gp = grid::gpar(fontsize = 15))
    dev.off()
  } else {
    pdf(paste(nameOut, "FinalModel.pdf", sep = "_"), width = dim(tg3)[2] + dim(tg1)[2] + dim(tg2)[2] + 14.5, height = dim(tg1)[1] / 2, useDingbats = F)
    gridExtra::grid.arrange(tg3, tg1, tg2, ncol = 3, nrow = 1, padding = 0, top = 0, left = 0)
    grid::grid.text(paste("Total variance explained: ", sum(as.numeric(varExpldepend)), "%", sep = ""), x = grid::unit(0.75, "npc"), y = grid::unit(0.90, "npc"), gp = grid::gpar(fontsize = 15))
    dev.off()
  }
  ### tb3out
  if (NetModelSel == "Adjusted") {
    nodeOut <- tb2out[, c(1, 4)]
  } else if (NetModelSel == "Univariate") {
    nodeOut <- tb3out[, c(1, 4)]
  } else {
    nodeOut <- tb2out[, c(1, 3)]
  }
  colnames(nodeOut)[2] <- "VarianceExplained"
  nodeOut$ID <- namesDf$newNames[match(nodeOut$Features, namesDf$originalNames)]
  nodeOut <- nodeOut[, c(3, 2, 1)]
  nodeOut$varexpl <- 1
  
  # other tgat are connected but not significant
  if (nrow(linkOut) > 0) {
    addAll <- data.frame(ID = unique(c(linkOut$from, linkOut$to)), VarianceExplained = 0)
    addAll$Features <- namesDf$originalNames[match(addAll$ID, namesDf$newNames)]
    
    addAll <- addAll[!addAll$ID %in% nodeOut$ID, ]
    if(nrow(addAll)>0){
      addAll$varexpl <- 2
    }
    #
    nodeOutAll <- rbind(nodeOut, addAll)
  } else {
    nodeOutAll <- nodeOut
  }
  # create igraph object
  net <- igraph::graph.data.frame(linkOut, nodeOutAll, directed = T)
  # rescale to size ans other attributes
  lsize <- rescale(igraph::V(net)$VarianceExplained, 0, 100, 0, 75)
  lsize[which(lsize > 0)] <- lsize[which(lsize > 0)] + 2
  igraph::V(net)$size <- lsize
  lcol <- rep("black", nrow(nodeOutAll))
  lcol[which(nodeOutAll$varexpl == 2)] <- "#B14D8E"
  igraph::V(net)$label.color <- lcol
  igraph::V(net)$label <- wrapNames(gsub(" 0%", "", paste(igraph::V(net)$Features, paste(igraph::V(net)$VarianceExplained, "%", sep = ""), sep = " ")), 8)
  
  ###introduce direction here to colour
  if(isTRUE(isads)){
    AnotaColours <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],RColorBrewer::brewer.pal(8,"Reds")[c(2,6)],RColorBrewer::brewer.pal(8,"Greens")[c(4,8)], RColorBrewer::brewer.pal(8,"Greens")[c(2,6)],RColorBrewer::brewer.pal(8,"Blues")[c(4,8)])
    names(AnotaColours) <- c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown","bufferingmRNAUp","bufferingmRNADown")
    coloursTmp <- as.character(AnotaColours[grep(regulationGen, names(AnotaColours))])
    #
    direct <- as.character()
    for(i in 1:nrow(nodeOutAll)){
      if(nodeOutAll[i,]$varexpl==1){
        IDtmp <- nodeOutAll[i,]$ID
        #
        directTmp <- as.numeric(cor.test(dataIn[,colnames(dataIn)==IDtmp,], dataIn$effM)[4])
        direct[i] <- ifelse(directTmp>0, 'plus','minus')
      } else {
        direct[i] <- 'nonsign'
      }
    }
    nodeOutAll$direct <- direct
    #
    colrs <- rep("white", nrow(nodeOutAll))
    colrs[which(nodeOutAll$direct == 'plus')] <- coloursTmp[1]
    colrs[which(nodeOutAll$direct == 'minus')] <- coloursTmp[2]
    igraph::V(net)$color <- colrs # colrs[igraph::V(net)]
  } else {
    colrs <- rep("#B0F2BC", nrow(nodeOutAll))
    colrs[which(nodeOutAll$varexpl == 2)] <- "white"
    igraph::V(net)$color <- colrs # colrs[igraph::V(net)]
  }
  # Set edges width based on weight:
  if (length(igraph::E(net)$weight) > 0) {
    if(isTRUE(useCorel)){
      igraph::E(net)$width <- rescale(abs(igraph::E(net)$weight), 0, 1, 0, 15)
    } else {
      igraph::E(net)$width <- rescale(abs(igraph::E(net)$weight), 0, 50, 0, 2.5)
    }
  }
      
  # change arrow size and edge color:
  igraph::E(net)$arrow.size <- .0
  ecolor <- rep(grDevices::adjustcolor("#F40009", alpha.f = 0.2), length(igraph::E(net)$width))
  ecolor[which(igraph::E(net)$weight < 0)] <- grDevices::adjustcolor("#1C39BB",alpha.f = 0.2)
  igraph::E(net)$color <- ecolor
        
  # We can also override the attributes explicitly in the plot:
  pdf(paste(nameOut, "network.pdf", sep = "_"), height = 8, width = 8, useDingbats = F,family = "Helvetica")
  par(bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 1.9, cex.lab = 1.5)
  par(mar = c(5, 5, 5, 5), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
  plot(net, shape = "sphere", vertex.label.font = 2, vertex.label.cex = 1, vertex.frame.color = "black", layout = layoutCalc(net, n = 2))
  if(isTRUE(isads)){
    tmpCoord <- legend("bottomleft", fill = coloursTmp, c("Positive regulation","Negative regulation"), cex = 1.1, bty = "n", xpd = T, inset = -0.1)
    legend(-1.5, 1.1, "Covary with significant features", bty="n",text.col="#B14D8E")
  } else {
    tmpCoord <-legend("bottomleft", fill = c("#B0F2BC"), "In Omnibus model", cex = 1.3, bty = "n", xpd = T, inset = -0.1)
  }
  legend(-1.5, 1.5, lwd = c(7, 3, 3, 7), col = c(rep(grDevices::adjustcolor("#F40009", alpha.f = 0.2), 2), rep(grDevices::adjustcolor("#1C39BB",alpha.f = 0.2), 2)), title = ifelse(isTRUE(useCorel),"Correlation", "Co-variance"), c("+", "", "", "-"), bty = "n", xpd = T)
  legend(-0.8, 1.5, pt.cex = c(4, 3, 2, 1), pch = 20, col = "gray75", title = c("Variance explained"), c("", "", "", ""), bty = "n", xpd = T)
  dev.off()
}
