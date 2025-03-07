anota2seqGetDirectedRegulations <- function(ads) {
  n_contrasts <- ncol(ads@contrasts)
  regModeList <- vector("list", length = n_contrasts)
  
  for (i in 1:n_contrasts) {
    translation_data <- ads@selectedTranslation@selectedRvmData[[i]]
    translated_data <- ads@selectedTranslatedmRNA@selectedRvmData[[i]]
    buffering_data <- ads@selectedBuffering@selectedRvmData[[i]]
    abundance_data <- ads@mRNAAbundance@totalmRNA[[i]]
    total_mrna_data <- ads@selectedTotalmRNA@selectedRvmData[[i]]
    
    translationUp <- translation_data[translation_data$apvEff > 0 & translation_data$singleRegMode == "translation", ]
    translationDown <- translation_data[translation_data$apvEff < 0 & translation_data$singleRegMode == "translation", ]
    translatedmRNAUp <- translated_data[translated_data$apvEff > 0, ]
    translatedmRNADown <- translated_data[translated_data$apvEff < 0, ]
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

regSel <- function(annot, region){#, ext=FALSE){
  nc <- grep(region, colnames(annot))
  #
  seqTmp <- annot[,nc]
  lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
  #
  #if(isTRUE(ext)){
  #  seq <- list()
  #  seq[[1]] <- annot$CDS_seq
  #  seq[[2]] <- annot$UTR3_seq
  #  extSeq <- combSeq(seqIn = seq)
  #  extSeq <- unlist(extSeq)
  #  #
  #  annotOut <- cbind(annot[,c(1:2)], seqTmp,lenTmp,extSeq)
  #} else {
  annotOut <- cbind(annot[,c(1:2)], seqTmp,lenTmp)
  #}
  #
  return(annotOut)
}

isoSel <- function(annot, method, setSeed){
  #Select per gene level
  if(method=='shortest'){
    annotOut <- as.data.frame(annot %>% group_by(geneID) %>% dplyr::slice(which.min(lenTmp)))
  } else if(method=='longest'){
    annotOut <- as.data.frame(annot %>% group_by(geneID) %>% dplyr::slice(which.max(lenTmp)))
  } else {
    if(!is.null(setSeed)){
      set.seed(setSeed)
    }
    annotOut <- as.data.frame(annot %>% group_by(geneID) %>% dplyr::slice_sample(n = 1))
  }
  return(annotOut)
}

adjustSeq <- function(annot,
                      adjObj,
                      region_adj,
                      excl = FALSE,
                      keepAll = FALSE) {
  
  checkAnnot(annot)
  #
  if(!is.null(adjObj)){
    check_adjObj(adjObj)
    valid_regions <- c('UTR5', 'UTR3')
    if (is.null(region_adj) | !all(region_adj %in% valid_regions)) {
      stop("'region_adj' has to be provided and can be only: 'UTR5','UTR3'. It should also match named entries in the list adjObj ")
    }
    if(!check_logical(excl)){
      stop("'excl' can only be only be logical: TRUE of FALSE ")
    }
    if(!check_logical(keepAll)){
      stop("'keepAll' can only be only be logical: TRUE of FALSE ")
    }
  }
  annotTmp <- annot
  #
  for(reg in region_adj){
    adjObj_temp <- adjObj[[reg]]

     if(is.null(adjObj_temp)){
      stop("The one or more regions specified in 'region_adj' do not match those provided in the 'adjObj' list. Please ensure the names of 'adjObj' are 'UTR5' and/or 'UTR3', and correspond to 'region_adj'.")
    }
 
    if(length(which(names(adjObj_temp) %in% annotTmp$id))==0){
      stop("It looks like transcript IDs provided in 'adjObj' do not match with transcript IDs in the existing annotation.")
    }
    adjObj_temp <- adjObj_temp[names(adjObj_temp) %in% annotTmp$id]
    #
    if (isTRUE(excl)) {
      annotTmp <- annotTmp[annotTmp$id %in% names(adjObj_temp), ]
    }
    #remove these that are not in annotation but in adjustment 
    #toRemove <- setdiff(names(adjObj_temp), annotTmp$id)
    #adjObj_temp <- adjObj_temp[!names(adjObj_temp) %in% toRemove]
    #
    if(length(adjObj_temp) == 0){
      stop("None of the entries in the adjustment vector are present in the existing annotation.")
    }
    #
    annotTmp[match(names(adjObj_temp), annotTmp$id), ifelse(reg=="UTR5","UTR5_seq","UTR3_seq")] <- adjObj_temp
    if (!isTRUE(keepAll)) {
      annotTmp <- annotTmp[(annotTmp$geneID %in% unique(annotTmp[annotTmp$id %in% names(adjObj_temp), ]$geneID) & annotTmp$id %in% unique(annotTmp[annotTmp$id %in% names(adjObj_temp), ]$id)) | (!annotTmp$geneID %in% unique(annotTmp[annotTmp$id %in% names(adjObj_temp), ]$geneID)), ]
    }
  }
  return(annotTmp)
}

resSel <- function(ads=NULL, 
                   regulation=NULL,
                   contrast, 
                   geneList){
  resOut <- list()
  #Extract all results
  if(!is.null(ads)){
    results <- anota2seqGetDirectedRegulations(ads)
    #
    if (!is.null(regulation)){
      res <- vector("list", length = length(regulation))
      if(!length(setdiff(contrast, seq(from = 1, to = dim(ads@contrasts)[2])))==0){
        stop('One or more of the contrasts provided are not included in the anota2seq object.')
      }
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
  } else {
    res <- geneList
  }
  resOut <- res[lapply(res,length)>3]
  return(resOut)
}

getBg <- function(ads=NULL, customBg=NULL, geneList=NULL){
  bgOut <- list()
  #Extract all results
  if(!is.null(ads)){
    bgOut <-  row.names(ads@dataP)
  } else {
    if(!is.null(customBg)){
      #add here to be sure that all genelist are in bg
      bgOut <- customBg
      if(!is.null(geneList)){
        tmpDiff <- setdiff(as.character(unlist(geneList)), bgOut)
        if(length(tmpDiff)>0){
          stop('There are genes in the geneList that are not in custom bg')
        }
      }
    } else {
      bgOut <- as.character(unlist(geneList))
    } 
  }
  return(bgOut)
}

check_id_type <- function(id) {
  # Check if the ID is purely numeric (likely an Entrez ID)
  if (grepl("^[0-9]+$", id)) {
    return("entrezID")
  }
  # Check if the ID is alphabetic or alphanumeric (likely a GeneID)
  else if (grepl("^[A-Za-z0-9]+$", id)) {
    return("geneID")
  } 
  # If it doesn't match either, return unknown
  else {
    return("unknown")
  }
}



coloursSel <- function(ads, genesIn, geneList, geneListcolours){
  coloursOut <- as.character()
  if(!is.null(ads)){
    AnotaColours <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],RColorBrewer::brewer.pal(8,"Reds")[c(2,6)],RColorBrewer::brewer.pal(8,"Greens")[c(4,8)], RColorBrewer::brewer.pal(8,"Greens")[c(2,6)],RColorBrewer::brewer.pal(8,"Blues")[c(4,8)])
    names(AnotaColours) <- c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown","bufferingmRNAUp","bufferingmRNADown")
    #
    coloursOut <- AnotaColours[gsub("\\_.*","",names(genesIn))]
  } else {
    coloursOut <-  geneListcolours
  } 
  return(coloursOut)
}

effectSel <- function(ads, regulationGen, contrastSel, effectMeasure){
  effM <- as.numeric()
  if (!is.null(effectMeasure)) {
    effM <- effectMeasure
  } else if (!is.null(ads)) {
    if(regulationGen=='mRNAAbundance'){
      regTmp <- 'totalmRNA'
    } else {
      regTmp <- regulationGen
    }
    scOut <- anota2seq::anota2seqGetOutput(ads, output = "singleDf", selContrast = contrastSel, getRVM = TRUE)
    effM <- scOut[, grepl(paste(regTmp, "apvEff", sep = "."), colnames(scOut))]
    names(effM) <- scOut$identifier
  } else {
    stop("Please provide anota2seq object or custom effect measure to be explained by provided features")
  }
  return(effM)
}

resQuant <- function(qvec, ptn){
  resOut <- list()
  if(!is.null(ptn_background(ptn))){
    res <- c(list(background=ptn_background(ptn)),ptn_geneList(ptn))
  } else {
    res <- ptn_geneList(ptn)
  }
  
  for(i in 1:length(res)){
    resOut[[names(res)[i]]] <- qvec[names(qvec) %in% res[[i]]]
  }
  return(resOut)
}

colPlot <- function(ptn){
  if(!is.null(ptn_background(ptn))){
    colOut <- c('grey45',ptn_colours(ptn))
  } else {
    colOut <- ptn_colours(ptn)
  }
  return(colOut)
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
      #yposTmp <- #range(as.numeric(unlist(resOut)))[2]#ifelse(range(as.numeric(unlist(resOut)))[2] <= 1, 1.1,range(as.numeric(unlist(resOut)))[2]+ j*1)
      
      rect(xleft = compTmp[1], xright = compTmp[2], ybottom = j-1  , ytop = j-1, lwd = 2)
      #
      #text(sum(compTmp) / 2, ifelse(range(as.numeric(unlist(resOut)))[2] <= 1,yposTmp + 0.05,yposTmp + 1), pvalTmp, cex = 0.75)
      text(sum(compTmp) / 2, j-0.5, pvalTmp, cex = 0.75)
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

plotUtils <- function(resOut, colOut, comparisons, ylabel, plotType) {
  
  if(plotType == "boxplot" | plotType == "violin"){
    m <- layout(mat = matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(1, 5))
    xlimTmp <- c(0.5, length(resOut) + 1.5)
    #
    par(mar = c(0, 8, 3, 0),bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    ylimTmp1 <- ifelse(!is.null(comparisons), length(comparisons), 0)
    plot(1,ylimTmp1, xlim=xlimTmp, ylim=c(0,ylimTmp1),xaxt = "n",type="n", yaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2, frame.plot = FALSE)
    if(!is.null(comparisons)){
      addStats(comparisons, plotType='boxplot', resOut, colOut)
    }
    #
    dataTmp <- as.numeric(unlist(resOut))
    ylimTmp2 <- roundNice(quantile(dataTmp,0.9), direction='up')
    par(mar = c(8, 8, 0, 0), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    plot(1,ylimTmp2, xlim=xlimTmp, ylim=c(0,ylimTmp2), xaxt = "n",type="n", yaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2, frame.plot = FALSE)
    #
    if(ylabel == 'Length (Log2 scale)'){
      axis(side = 2, font = 2, las = 2, lwd = 2, at = sapply(c(1, 25, 100, 200, 400, 1000, 4000, 25000), log2), labels = c(0, 25, 100, 200, 400, 1000, 4000, 25000))
      mtext(side = 2, line = 6,  ylabel, col = "black", font = 2, cex = 1.7, at = median(dataTmp))
    } else {
      axis(side = 2, font = 2, las = 2, lwd = 2)
      mtext(side = 2, line = 6,  ylabel, col = "black", font = 2, cex = 1.7, at = roundNice(median(dataTmp), direction='up'))
    }
    text(1:length(resOut), par("usr")[3] - 0.45, labels = names(resOut), xpd = NA, cex = 0.9, srt = 45, adj = 1)
    #
    if (names(resOut)[1] == 'background') {
      abline(lty = 5, h = median(resOut[[1]]))
    }
    #
    for (i in 1:length(resOut)) {
      if(plotType=='boxplot'){
        boxplot(resOut[[i]], add = TRUE, at = i, col = colOut[i], xaxt = "n", xlab = "", ylab = "", type = "n", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE, outcol = "grey65", whiskcol = "grey65", outline = FALSE, medcol = "black", staplelty = 0, whisklty = 1)
      } else if (plotType=='violin'){
        vioplot::vioplot(resOut[[i]], add = TRUE, at = i, col = colOut[i], xaxt = "n", xlab = "", ylab = "", main = "", lwd = 1, bty = "n", yaxt = "n", font = 2, frame.plot = FALSE)
      } 
      if(ylabel == 'Length (Log2 scale)'){
        text(i, 0, round(mean(antilog(resOut[[i]], 2), 0)), font = 2)
      } else {
        text(i, 0, round(mean(resOut[[i]]),0), font = 2)
      }
    } 
  } else if(plotType == "ecdf"){
    #xlim_min <- roundNice(as.numeric(quantile(as.numeric(unlist(resOut)), 0.01)),direction='down')
    #xlim_max <- roundNice(as.numeric(quantile(as.numeric(unlist(resOut)), 0.99)),direction='up')
    
    par(mar = c(5, 5, 8, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
    plot(ecdf(resOut[[1]]), col = colOut[1], main = "", xlab = "", ylab = "", verticals = TRUE, do.p = FALSE, lwd = 3, bty = "n", font = 2)#, yaxt = "none", xlim = c(xlim_min, xlim_max), xaxt = "none")
    
    mtext(side = 1, line = 4, ylabel, col = "black", font = 2, cex = 1.2)
    mtext(side = 2, line = 3, "Fn(x)", col = "black", font = 2, cex = 1.2)
    
    #axis(side = 1, seq(xlim_min, xlim_max, 1), font = 2, lwd = 2)
    #axis(side = 2, seq(0, 1, 0.2), font = 2, las = 2, lwd = 2)
    for (i in 2:length(resOut)) {
      lines(ecdf(resOut[[i]]), col = colOut[i], main = "", xlab = "", verticals = TRUE, do.p = FALSE, lwd = 4)
    }
    if(!is.null(comparisons)){
      addStats(comparisons, plotType='ecdf', resOut, colOut)
    }
  }  else {
      stop("Please provide correct 'plotType'")
  }
}

antilog<-function(lx,base){ 
  lbx<-lx/log(exp(1),base=base) 
  result<-exp(lbx) 
  result 
} 

is_by_3 <- function(seqs) {
  all(sapply(seqs, function(x) length(seqinr::c2s(x)) %% 3 == 0))
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
    #gROut <- GenomicRanges::reduce(GenomicRanges::GRanges(seqnames='tmp', ranges=IRanges::IRanges(start=motOut,end=motOut)),min.gapwidth=len)
    #
    dtTmp <- data.table::data.table(start = motOut, end = motOut)
    data.table::setorder(dtTmp, start)
    dtTmp[, group := cumsum(c(1, diff(start) > len))]
    dtTmp <- dtTmp[, .(start = min(start), end = max(end)), by = group]

    if(unit == 'number'){
      #nMot <- length(gROut)
      nMot <- nrow(dtTmp)
    } else if (unit == 'position') {
      nMot <- list()
      #nMot[["start"]] <- as.numeric(start(gROut@ranges))
      #nMot[["end"]] <- as.numeric(end(gROut@ranges)) + lenTmp - 1
      nMot[["start"]] <- as.numeric(dtTmp$start)
      nMot[["end"]] <- as.numeric(dtTmp$end)  + lenTmp - 1
    }
  } else {
    nMot <- ifelse(unit == "number", 0, NA)
  }
  return(nMot)
}

#Combine 
calc_g4 <- function(x, min_score, unit){
    seqTmp <- DNAString(x)
    predTmp <- pqsfinder::pqsfinder(seqTmp, min_score = min_score, strand = '+')
    if(nrow(predTmp@elementMetadata)>0){
      if(unit == 'number'){
        nMot <- nrow(predTmp@elementMetadata)
      } else if (unit == 'position') {
        nMot <- list()
        #nMot[["start"]] <- as.numeric(start(gROut@ranges))
        #nMot[["end"]] <- as.numeric(end(gROut@ranges)) + lenTmp - 1
        nMot[["start"]] <- as.numeric(predTmp@ranges@start)
        nMot[["end"]] <- as.numeric(predTmp@ranges@start) + as.numeric(predTmp@ranges@width) -1
      }
    } else {
      nMot <- ifelse(unit == "number", 0, NA)
    }
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

remove_last3 <- function(seq) {
  seqOut <- sapply(seq, function(x) {
    #
    substring(x, 1, nchar(x) - 3)
  })
  return(seqOut)
}


#
codonCount <- function(seq, gene, codonN=1){
  #
  if(codonN==1){
    #
    tmpEff <- seqinr::uco(seqinr::s2c(seq),index = "eff")
    tmpFreq <- seqinr::uco(seqinr::s2c(seq),index = "freq")
    #
    tmpCodon <- data.frame(geneID=gene,codon=toupper(names(tmpEff)),AA=seqinr::translate(seqinr::s2c(seqinr::c2s(toupper(names(tmpEff))))), count=as.numeric(tmpEff),frequency=as.numeric(tmpFreq))
    tmpCodon <- tmpCodon %>% mutate(AA = case_when(codon == "TGA" ~ "U",codon == "TAG" ~ "O",codon == "TAA" ~ "Stop", TRUE ~ AA))
    
    tmpCodon <- tmpCodon %>% group_by(AA) %>% mutate(AACountPerGene=sum(count))
    tmpCodon$relative_frequency <- tmpCodon$count/tmpCodon$AACountPerGene
    tmpCodon$relative_frequency[is.na(tmpCodon$relative_frequency)] <- 0
  } else if(codonN > 1){
    seqIn <- seqinr::s2c(tolower(seq))
    
    tmpEff <- seqinr::count(seq = seqIn, wordsize = 3*codonN, start = 3, by = 3,freq = FALSE)
    tmpFreq <- tmpEff/sum(tmpEff)
    
    tmpCodon <- data.frame(geneID=gene,codon=toupper(names(tmpEff)), count=as.numeric(tmpEff), AA=NA, frequency=as.numeric(tmpFreq),AACountPerGene=NA, relative_frequency=NA)
  }
  #
  return(tmpCodon)
}


#divide_3 <- function(stringIn) {
#  # Check if the string length is divisible by 3
#  if (nchar(stringIn) %% 3 != 0) {
#    stop("The length of the string is not divisible by 3.")
#  }
#  #
#  tmpChar <- strsplit(stringIn, "")[[1]]
#  #
#  tmpLen <- nchar(stringIn)
#  
##  # Divide the string into groups of three
#  out <- sapply(seq(1, tmpLen, by = 3), function(i) {
#    paste(tmpChar[i:min(i+2, tmpLen)], collapse = "")
#  })
#  
#  return(out)
#}

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


#roundNice <- function(x, nice=c(1,2,4,5,6,8,10), direction) {
#     if(length(x) != 1) stop("'x' must be of length 1")
#     if(direction == 'up'){
#      10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
#     } else if (direction == 'down'){
#       10^floor(log10(x))
#     } else {
#       stop('wrong input for rounding function')
#     }
#}

roundNice <- function(x, nice=c(1,2,4,5,6,8,10), direction) {
  if(length(x) != 1) stop("'x' must be of length 1")
  
  # Handling the sign of the input
  sign_x <- sign(x)
  x_abs <- abs(x)
  
  if(direction == 'up'){
    result <- 10^floor(log10(x_abs)) * nice[[which(x_abs <= 10^floor(log10(x_abs)) * nice)[[1]]]]
  } else if (direction == 'down'){
    result <- 10^floor(log10(x_abs)) * nice[[tail(which(x_abs >= 10^floor(log10(x_abs)) * nice), 1)]]
  } else {
    stop('wrong input for rounding function')
  }
  
  return(sign_x * result)
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

dup_tolerance <- function(numIn, tol = 1e-8) {
  sapply(seq_along(numIn), function(i) {
    any(abs(numIn[i] - numIn[-i]) < tol)
  })
}

getDup <- function(numIn, tol = 1e-8) {
  #
  dupsOut <- list()

  #
  checked <- rep(FALSE, length(numIn))
  
  #
  for (i in seq_along(numIn)) {
    if (!checked[i]) {
      #
      dups <- which(abs(numIn - numIn[i]) < tol)
      
      #
      if (length(dups) > 1) {
        dupsOut[[length(dupsOut) + 1]] <- dups
      }
      checked[dups] <- TRUE
    }
  }
  return(dupsOut)
}


printDup <- function(dupIn) {
  for (i in seq_along(dupIn)) {
    namesMerged <- paste(names(dupIn[[i]]), collapse = " and ")
    messageOut <- paste("These elements are too corelated to each other: ", namesMerged)
    print(messageOut)
  }
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

s4_to_dataframe <- function(obj) {
  tmpNames <- slotNames(obj)
  
  tmpOut <- lapply(tmpNames, function(x) slot(obj, x))
  out <- as.data.frame(setNames(tmpOut, tmpNames))
  return(out)
}


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

prepFeatures <- function(ptn, 
                         features){
  check_ptn(ptn)
  if (!is_valid_named_list(features)){
    stop("features should be a named list of numeric vectors only")
  }
  #effM <- ptn_eff(ptn)
  #if(!is_numeric_vector(ptn_eff(ptn))){
  #  stop("'effectMeasure' should be a numeric vector")
  #}
  #
  featureNames <- names(features)
  #featuresTmp <- append(features, list(effM))
  #featuresTmp <- unname(featuresTmp)
  
  #tmpDf <- data.frame(t(plyr::ldply(featuresTmp, rbind,.id = NULL)))
  tmpDf <- data.frame(t(plyr::ldply(features, rbind,.id = NULL)))
  #colnames(tmpDf) <- c(featureNames,'effM')
  colnames(tmpDf) <- featureNames
  
  datOut <- na.omit(tmpDf)
  message(paste(nrow(tmpDf)-nrow(datOut), 'genes removed because of NAs', sep=' '))
  
  ptn@features <- datOut
  return(ptn)
}

normalizeLayout <- function(layout) {
  layout[, 1] <- (layout[, 1] - min(layout[, 1])) / (max(layout[, 1]) - min(layout[, 1])) * 2 - 1
  layout[, 2] <- (layout[, 2] - min(layout[, 2])) / (max(layout[, 2]) - min(layout[, 2])) * 2 - 1
  return(layout)
}


runLM <- function(dataIn, namesDf, allFeat, useCorel, covarFilt, nameOut, NetModelSel, coloursIn, lmfeatGroup, lmfeatGroupColour){
  #
  fval <- list()
  pval <- list()
  
  #Univariate
  models <- lapply(colnames(dataIn)[-length(colnames(dataIn))], function(x) {
    anova(lm(substitute(effM ~ i, list(i = as.name(x))), data = dataIn))
  })
  names(models) <- namesDf$originalNames
  #
  step1expl <- round(sapply(models, function(x) (x[1, 2] / (sum(x[1, 2], x[2, 2]))) * 100), 2)
  step1pval <- sapply(models, function(x) (x[1, 5]))
  step1fdr <- p.adjust(step1pval)
  
  #
  uniOut <- new("postNetUnivariate",
                pvalue = step1pval,
                fdr = step1fdr,
                varianceExplained = step1expl
  )
  
  presel <- as.numeric(which(step1fdr > 0.05 | is.na(step1fdr)))

  if (length(presel) > 0) {
    step1expl <- step1expl[-presel]
    step1pval <- step1pval[-presel]
    step1fdr <- step1fdr[-presel]
  }

  if (!isTRUE(allFeat) & length(presel) > 0) {
    dataIn <- dataIn[, -presel]
    models <- models[-presel]
  }

  stepwiseModels <- list()
  #
  fvalTmp <- sapply(models, function(x) x[1, 4])
  names(fvalTmp) <- colnames(dataIn)[-length(colnames(dataIn))]
  #
  if(length(getDup(fvalTmp))>0){
    names(fvalTmp) <- getOrgNames(names(fvalTmp),namesDf)
    dupTmp <- getDup(fvalTmp)
    printDup(dupTmp)
    stop('Two of the input features are almost perfectly correlated. Only one of the correlated features should be included')
  }
  fval[[1]] <- fvalTmp

  pvalTmp <- sapply(models, function(x) x[1, 5])
  names(pvalTmp) <- colnames(dataIn)[-length(colnames(dataIn))]
  pval[[1]] <- pvalTmp

  bestSel <- names(which.max(fvalTmp))
  outSel <- names(which(pvalTmp > 0.05 | is.na(pvalTmp)))

  i <- 1
  # 
  while (length(colnames(dataIn)[-length(colnames(dataIn))][!colnames(dataIn)[-length(colnames(dataIn))] %in% c(bestSel, outSel)]) > 0) {
    #
    i <- i + 1
    #
    models2 <- list()
    #
    x_sel <- paste(bestSel, collapse = " + ")
    #
    for (j in 1:(length(colnames(dataIn)[-length(colnames(dataIn))]) - length(c(bestSel, outSel)))) {
      varx <- colnames(dataIn)[-length(colnames(dataIn))][!colnames(dataIn)[-length(colnames(dataIn))] %in% c(bestSel, outSel)][j]
      design <- as.formula(paste(paste("effM", paste(x_sel, collapse = " + "), sep = "~"), varx, sep = "+"))
      models2[[j]] <- anova(lm(design, data = dataIn))

      isNAout <- names(lm(design, data = dataIn)[[1]][which(is.na(lm(design, data = dataIn)[[1]]))])
      
      if(length(isNAout)>0){
        stop(paste('Consider removing: ',getOrgNames(names(isNAout),namesDf),' as it is too correlated with other features'))
      }
      rownames(models2[[j]]) <- c(getOrgNames(rownames(models2[[j]])[-length(rownames(models2[[j]]))],namesDf),"Residuals")
    }
    #
    stepwiseModels[[paste('step', i, sep=' ')]] <- models2
    #
    fvalTmp <- sapply(models2, function(x) x[nrow(x) - 1, 4])
    names(fvalTmp) <- colnames(dataIn)[-length(colnames(dataIn))][!colnames(dataIn)[-length(colnames(dataIn))] %in% c(bestSel, outSel)]
    #
    if(length(getDup(fvalTmp))>0){
      names(fvalTmp) <- getOrgNames(names(fvalTmp), namesDf)
      dupTmp <- getDup(fvalTmp)
      printDup(dupTmp)
      stop('Two of the input features are almost perfectly correlated. Only one of the correlated features should be included')
    }
    #
    pvalTmp <- sapply(models2, function(x) x[nrow(x) - 1, 5])
    names(pvalTmp) <- colnames(dataIn)[-length(colnames(dataIn))][!colnames(dataIn)[-length(colnames(dataIn))] %in% c(bestSel, outSel)]
    #

    fval[[i]] <- fvalTmp
    pval[[i]] <- pvalTmp

    bestTmp <- names(which.max(fvalTmp))
    outTmp <- names(which(pvalTmp > 0.05 | is.na(pvalTmp)))
    #
    if (length(outTmp) > 0) {
      outSel <- append(outSel, outTmp)
      if( length(bestTmp) > 0) {
        if (!bestTmp %in% outTmp) {
          bestSel <- append(bestSel, bestTmp)
        }
      }
    } else {
      if( length(bestTmp) > 0) {
        bestSel <- append(bestSel, bestTmp)
      }
    }
  }
  #
  tmp_fval <- t(plyr::ldply(fval, rbind))
  tmp_fval <- tmp_fval[c(bestSel, rev(outSel)), ]

  tmp_pval <- t(plyr::ldply(pval, rbind))
  tmp_pval <- tmp_pval[c(bestSel, rev(outSel)), ]
  
  tb1Out <- round(apply(tmp_fval, 2, as.numeric), 2)
  row.names(tb1Out) <- namesDf$originalNames[match(row.names(tmp_fval), namesDf$newNames)]
  tb1Out[is.na(tb1Out)] <- ""
  colnames(tb1Out) <- paste("step", seq(1, ncol(tb1Out), 1), sep = "")

  colours <- matrix("white", nrow(tb1Out), ncol(tb1Out))
  linkIn <- matrix(NA, nrow(tb1Out), ncol(tb1Out))
  
  for (k in 2:nrow(colours)) {
    trow <- as.numeric(na.omit(as.numeric(tb1Out[k, ])))
    #
    if(length(trow)>0){
      tsel_tab <- which((lag(trow) / trow) > 2)
      tsel_net <- lag(trow) - trow
      #
      colours[k, tsel_tab] <- "#FDE0C5"
      linkIn[k, 1:length(tsel_net)] <- tsel_net
    }
  }
  
  for (k in 1:ncol(colours)) {
    colours[k, k] <- "#B0F2BC"
  }
  #
  colours[which(tmp_pval > 0.05)] <- "#B14D8E"
  #colours <- colours[apply(colours, 1, function(x) any(x != "white")), ]
  linkIn[which(abs(linkIn) < covarFilt)] <- NA
  row.names(linkIn) <- row.names(tmp_fval)
  #
  tt1 <- gridExtra::ttheme_default(core = list(fg_params = list(fontface = c(rep("plain", ncol(tb1Out)))), bg_params = list(fill = colours, col = "black")))#,colhead=list(fg_params=list(rot=90,hjust=0, y=0)))
  tg1 <- gridExtra::tableGrob(tb1Out, theme = tt1)
  #
  stepWiseout <- new("postNetStepWise",
                     models = stepwiseModels,
                     table = tb1Out)
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
    tmpM <- anova(lm(design, data = dataIn))

    varExplIndepend[i] <- round((tmpM[nrow(tmpM) - 1, 2] / sum(tmpM$`Sum Sq`)) * 100, 2)
  }
  names(varExplIndepend) <- bestSel
  
  #
  #varExplIndepend2 <- numeric()
  #
  #if (length(presel) > 0) {
  #  step1sel <- namesDf$newNames[-presel]
  #} else {
  #  step1sel <- namesDf$newNames
  #}
  #for (i in 1:length(step1sel)) {
  #  tmpFeat2 <- step1sel[i]
  #  restFeat2 <- step1sel[-i]
  #  
  #  design2 <- as.formula(paste(paste("effM", paste(restFeat2, collapse = " + "), sep = "~"), tmpFeat2, sep = "+"))
  #  tmpM2 <- anova(lm(design2, data = dataIn))
  #
  #  varExplIndepend2[i] <- round((tmpM2[nrow(tmpM2) - 1, 2] / sum(tmpM2$`Sum Sq`)) * 100, 2)
  #}
  #names(varExplIndepend2) <- step1sel
  
  #
  tb2out <- data.frame(Features = namesDf$originalNames[match(bestSel, namesDf$newNames)], Pvalue = format(varExpl$`Pr(>F)`[1:length(bestSel)], scientific = T, digits = 2), VarianceExplained_Omnibus = as.numeric(varExpldepend), VarianceExplained_Adjusted = as.numeric(varExplIndepend))
  tg2 <- gridExtra::tableGrob(tb2out, rows = NULL)
  
  rownames(tmpM) <- c(namesDf$originalNames[match(rownames(tmpM)[-length(rownames(tmpM))], namesDf$newNames)],"Residuals")
  
  finalModelout <- new("postNetFinalModel",
                       totalVarianceExplained = sum(as.numeric(varExpldepend)),
                       finalModel = tmpM,
                       table = tb2out)

  
  tb3out <- data.frame(Features = names(step1expl), Pvalue_Univariate = format(as.numeric(step1pval), scientific = T, digits = 2), FDRvalue_Univariate = format(as.numeric(step1fdr), scientific = T, digits = 2), VarianceExplained_Univariate = as.numeric(step1expl))
  tb3out <- tb3out[with(tb3out, order(-tb3out$VarianceExplained_Univariate)), ]
  
  tg3 <- gridExtra::tableGrob(tb3out, rows = NULL)
  
  # 
  #tb4out <- data.frame(Features = namesDf$originalNames[match(names(varExplIndepend2), namesDf$newNames)], VarianceExplained_IndependentAll = as.numeric(varExplIndepend2))
  #tb4out <- tb4out[with(tb4out, order(-tb4out$VarianceExplained_IndependentAll)), ]
  #tg4 <- gridExtra::tableGrob(tb4out, rows = NULL)
  
  #pdf(paste(nameOut, "varexpl_independAll.pdf", sep = "_") , width = dim(tg4)[2] + dim(tg4)[2] + 4, height = dim(tg4)[1] / 2, useDingbats = F)
  #gridExtra::grid.arrange(tg4, ncol = 1, nrow = 1)
  #dev.off()
  
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
      f1dat <- dataIn[,colnames(dataIn)==f1tmp]
      f2dat <- dataIn[,colnames(dataIn)==f2tmp]
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
  ###
  if (NetModelSel == "adjusted") {
    nodeOut <- tb2out[, c(1, 4)]
  } else if (NetModelSel == "univariate") {
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
  net <- igraph::graph.data.frame(linkOut, nodeOutAll, directed = F)
  # rescale to size ans other attributes
  lsize <- rescale(igraph::V(net)$VarianceExplained, 0, 100, 0, 75)
  lsize[which(lsize > 0)] <- lsize[which(lsize > 0)] + 2
  igraph::V(net)$size <- lsize
  lcol <- rep("black", nrow(nodeOutAll))
  lcol[which(nodeOutAll$varexpl == 2)] <- "#B14D8E"
  igraph::V(net)$label.color <- lcol
  igraph::V(net)$label <- wrapNames(gsub(" 0%", "", paste(igraph::V(net)$Features, paste(igraph::V(net)$VarianceExplained, "%", sep = ""), sep = " ")), 8)
  
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
  colrs[which(nodeOutAll$direct == 'plus')] <- coloursIn[1]
  colrs[which(nodeOutAll$direct == 'minus')] <- coloursIn[2]
  igraph::V(net)$color <- colrs # colrs[igraph::V(net)]

  if (length(igraph::E(net)$weight) > 0) {
    if(isTRUE(useCorel)){
      igraph::E(net)$width <- rescale(abs(igraph::E(net)$weight), 0, 1, 0, 15)
    } else {
      igraph::E(net)$width <- rescale(abs(igraph::E(net)$weight), 0, 50, 0, 2.5)
    }
  }
      
  # 
  igraph::E(net)$arrow.size <- .0
  ecolor <- rep(grDevices::adjustcolor("#F40009", alpha.f = 0.2), length(igraph::E(net)$width))
  ecolor[which(igraph::E(net)$weight < 0)] <- grDevices::adjustcolor("#1C39BB",alpha.f = 0.2)
  igraph::E(net)$color <- ecolor
  
  if(!is.null(lmfeatGroup)){
    igraph::E(net)$weight <- 1
    
    igraph::V(net)$Group <- as.vector(lmfeatGroup[match(igraph::V(net)$name, names(lmfeatGroup))])
  #
    for(i in unique(igraph::V(net)$Group)) {
      GroupV <- which(igraph::V(net)$Group == i)
      net <- igraph::add_edges(net, combn(GroupV, 2), attr=list(weight=5))
    } 
    layOut <-  igraph::layout.fruchterman.reingold(net,weights=igraph::E(net)$weight)
  }
  layOut <- layoutCalc(net, n = 2)
  layOut<- normalizeLayout(layOut)
  
  #
  pdf(paste(nameOut, "network.pdf", sep = "_"), height = 8, width = 8, useDingbats = F,family = "Helvetica")
  m <- layout(mat = matrix(c(1, 2, 3), nrow = 3, ncol = 1), heights = c(2, 8, 1))
  #
  par(mar = c(0, 5, 5, 5),bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)

  plot(-2,2,xlim=c(-2,2), ylim=c(-2,2), xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2, frame.plot = FALSE,xaxt = "n",type="n", yaxt = "n")
  legend(-2, 2, lwd = c(7, 3, 3, 7), col = c(rep(grDevices::adjustcolor("#F40009", alpha.f = 0.2), 2), rep(grDevices::adjustcolor("#1C39BB",alpha.f = 0.2), 2)), title = ifelse(isTRUE(useCorel),"Correlation", "Co-variance"), c("+", "", "", "-"), bty = "n", xpd = T)
  legend(-1.5, 2, pt.cex = c(4, 3, 2, 1), pch = 20, col = "gray75", title = c("Variance explained"), c("", "", "", ""), bty = "n", xpd = T)
  legend(0.5, 2, paste('Total variance explained: ', sum(igraph::V(net)$VarianceExplained), '%'),cex=1.25, col='grey30')
  
  
  par(bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 1.9, cex.lab = 1.5)
  par(mar = c(2, 5, 1, 5), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
  plot(net, vertex.label.font = 2, vertex.label.cex = 1, vertex.frame.color = "black", edge.curved = FALSE,rescale = FALSE,layout = layOut)#,shape = "sphere")#,layoutCalc(net, n = 2))
  
  if(!is.null(lmfeatGroup)){
    for (group in unique(igraph::V(net)$Group)) {
      #
      group_nodes <- which(igraph::V(net)$Group == group)
      # 
      group_coords <- layOut[group_nodes, ]
      #
      if (is.vector(group_coords)) {
        group_coords <- matrix(group_coords, nrow = 1)
      }
      #
      group_center <- colMeans(group_coords)
    
      #
      layout_range_x <- max(layOut[, 1]) - min(layOut[, 1])
      layout_range_y <- max(layOut[, 2]) - min(layOut[, 2])
      margin_x <- 0.075 * layout_range_x
      margin_y <- 0.075 * layout_range_y
      x_radius <- (max(group_coords[, 1]) - min(group_coords[, 1])) / 2 + margin_x
      y_radius <- (max(group_coords[, 2]) - min(group_coords[, 2])) / 2 + margin_y

      #
      plotrix::draw.ellipse(x = group_center[1], y = group_center[2], a = x_radius, b = y_radius,border = adjustcolor(lmfeatGroupColour[group],alpha.f = 0.5), lwd = 10)
      
      text_x <- group_center[1] + x_radius * 0.75
      text_y <- group_center[2] + y_radius * 0.75
      
      # Add group label at the calculated text position
      text(x = text_x, y = text_y, labels = paste(sum(igraph::V(net)$VarianceExplained[which(igraph::V(net)$Group == group)]),'%',sep=''), col = lmfeatGroupColour[names(lmfeatGroupColour) == group], cex = 1.25, font = 2, pos = 4)
      legend('topright', bty='n', names(lmfeatGroupColour), title='Groups',text.col = lmfeatGroupColour, title.col = 'grey30')
    }

    #varexplTmp <- as.numeric()
    #for (i in 1:length(unique(igraph::V(net)$Group))) {
    #  varexplTmp[i] <- sum(igraph::V(net)$VarianceExplained[which(igraph::V(net)$Group == unique(igraph::V(net)$Group)[i])])
    #}
    #names(varexplTmp)<- unique(igraph::V(net)$Group)
    #legend(1, 1, paste(varexplTmp, '%',sep=''), cex=1.25, text.col=lmfeatGroupColour, bty = 'n')
  }
  
  par(mar = c(5, 5, 0, 5),bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.4, cex.main = 1.7, cex.lab = 1.3)
  plot(-2,2,xlim=c(-2,2), ylim=c(-2,2), xlab = "", ylab = "", main = "", lwd = 1, bty = "n", font = 2, frame.plot = FALSE,xaxt = "n",type="n", yaxt = "n")
  tmpCoord <- legend(-2,-2, fill = coloursIn, c("Positive regulation","Negative regulation"), cex = 1.1, bty = "n", xpd = T, inset = -0.1)
  legend('topright', "Covary with significant features", bty="n",text.col="#B14D8E")
  dev.off()
  
  selectedFeatures <- nodeOutAll[nodeOutAll$varexpl==1,]$VarianceExplained
  names(selectedFeatures) <- nodeOutAll[nodeOutAll$varexpl==1,]$Features
  
  lmOut <- new("postNetFeatureIntegration_lm",
               univariateModel = uniOut,
               stepwiseModel = stepWiseout,
               finalModel = finalModelout,
               selectedFeatures = selectedFeatures,
               networkGraph = net
  )
  return(lmOut)
}




plotScatterInd <- function(set1, set2=NULL, orgName, coloursIn, nameOut ){
  
  pdf(paste(nameOut, orgName, "individually.pdf", sep = "_"), width = 8, height = 8, useDingbats = F)
  par(mar = c(9, 5, 5, 4), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.7, cex.main = 1.7, cex.lab = 1.3)
  #
  if(is.null(set2)){
    set <- set1
  } else {
    set <- rbind(set1, set2)
  }
  
  xlim_min <- roundNice(min(set[,1]),direction = 'down')
  xlim_max <- roundNice(max(set[,1]),direction = 'up')
  ylim_min <- roundNice(min(set[,2]),direction = 'down')
  ylim_max <- roundNice(max(set[,2]),direction = 'up')
  
  plot(set1[,1], set1[,2], col = coloursIn[1], xlim=c(xlim_min,xlim_max), ylim=c(ylim_min,ylim_max), pch = 16, cex = 1, xlab = "", ylab = "", lwd = 1, bty = "n", font = 2)
  if(!is.null(set2)){
    points(set2[,1], set2[,2], pch = 16, col = coloursIn[2])
  }
  #
  mtext(side = 2, line = 3, 'effM', col = "black", font = 2, cex = 1.7)
  #axis(side = 2, seq(-ylim, ylim, 2), font = 2, las = 2, lwd = 2)

  mtext(side = 1, line = 4, orgName, col = "black", font = 2, cex = 1.7, at = (xlim_min + xlim_max) / 2)
  #axis(side = 1, seq(xlim_min, xlim_max, ifelse((xlim_max - xlim_min) / 5 >= 0, roundUpNice((xlim_max - xlim_min) / 5), -roundUpNice(abs((xlim_max - xlim_min) / 5)))), font = 2, lwd = 2)
  #
  if (length(unique(set[,1])) > 2 & IQR(set[,1]) > 0) {
    f1 <- predict(smooth.spline(set[,2] ~ set[,1]))
    lines(f1$x[which(f1$x > xlim_min & f1$x < xlim_max)], f1$y[which(f1$x > xlim_min & f1$x < xlim_max)], col = "#AFBADC", lwd = 4, lend = 2)
    lines(f1$x[which(f1$x > xlim_min & f1$x < xlim_max)], f1$y[which(f1$x > xlim_min & f1$x < xlim_max)], col = "black", lwd = 1, lend = 2, lty = 3)
  }
  #
  if (!is.na(as.numeric(coefficients(lm(set[,2] ~ set[,1]))[2]))) {
    plotrix::ablineclip(lm(set[,2] ~ set[,1]), col = "#AFBADC", lwd = 4, x1 = xlim_min, x2 = xlim_max)
    plotrix::ablineclip(lm(set[,2] ~ set[,1]), col = "black", lwd = 1, x1 = xlim_min, x2 = xlim_max)
  
    text((xlim_min + xlim_max) / 2, ylim_max, paste("pvalue ", format(as.numeric(cor.test(set[,1], set[,2])[3]), scientific = T, digits = 3), ", r=", round(as.numeric(cor.test(set[,1], set[,2])[4]), 3), sep = ""), bty = "n", col = "black", cex = 1.25, font = 2)
  }
  dev.off()
}

getOrgNames <- function(newN, namesDf){
    outN <- namesDf$originalNames[match(newN, namesDf$newNames)]
    return(outN)
}

is_binary <- function(x) {
  if (!is.numeric(x)) {
    return(FALSE)
  }
  unique_values <- unique(x)
  if (length(unique_values) == 2 && all(unique_values %in% c(0, 1))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

plot_fmap <- function(fMap, colVec, remExtreme = NULL, name){
  
  if(!is.null(remExtreme)){
    minV <- quantile(colVec, remExtreme)
    maxV <- quantile(colVec, 1 - remExtreme)
    
    fmapRes$colVecColour <- pmin(pmax(colVec, minV), maxV)
  } else {
    fmapRes$colVecColour <- colVec
  }
  if(!is_binary(colVec)){
    colVecPlot <- ggplot2::ggplot(fmapRes, ggplot2::aes(x = fUMAP1, y = fUMAP2, color = colVecColour)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_gradient2(low = "blue",mid = 'white', high = "red") +
      ggplot2::labs(title = paste('Feature: ', name, sep=''),  x = "fUMAP 1", y = "fUMAP 2", color = name) +
      ggplot2::theme_minimal() +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") 
  
    #keep original legend
    tmpLeg <- data.frame(x = 1, y = 1,colVec = seq(min(colVec), max(colVec), length.out = 100))
  
    colVecLeg<- ggplot2::ggplot(tmpLeg, ggplot2::aes(x = x, y = y, color = colVec)) +
      ggplot2::geom_point(size = 0) +  # Invisible points
      ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red") +
      ggplot2::labs(color = name) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.key.height = ggplot2::unit(1.5, "cm"),legend.key.width = ggplot2::unit(0.75, "cm"))
  
    legend_grob <- ggplot2::ggplotGrob(colVecLeg)
    legendOut <- legend_grob$grobs[[which(legend_grob$layout$name == "guide-box")]]
  } else {
    colVecPlot <- ggplot2::ggplot(fmapRes, ggplot2::aes(x = fUMAP1, y = fUMAP2, color = factor(colVecColour))) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(values = c("0" = "grey75", "1" = "firebrick1")) +  # Custom color for binary values
      ggplot2::labs(title = paste('Feature: ', name, sep=''),  x = "fUMAP 1", y = "fUMAP 2", color = name) +
      ggplot2::theme_minimal() +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
    
    
    legDataTmp<- data.frame(category = name)
    colVecLeg <- ggplot2::ggplot(legDataTmp, ggplot2::aes(x = 1, y = 1)) +
      ggplot2::geom_point(size = 4, color = "firebrick1") +
      ggplot2::labs(color = name) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold"))
    
    legend_grob <-  ggplot2::ggplotGrob(colVecLeg)
    legendOut <- legend_grob
    
  }
  
  return(list(mainPlot = colVecPlot, legend = legendOut))  
}
