goAnalysis <- function(a2sU,
                       #ads=NULL,
                       #regulation=NULL,
                       #contrast=NULL,
                       genesSlopeFiltOut = NULL,
                       #geneList = NULL,
                       #customBg = NULL,
                       #species,
                       category, #To choose from (gene ontologies) 'BP', 'CC', 'MF' or 'KEGG'
                       maxSize = 500,
                       minSize = 10,
                       counts = 10,
                       FDR = 0.15,
                       name = NULL
){
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  species <- a2sU_species(a2sU)
  if (!species %in% c("human","mouse")) {
    stop("This option is only  available for species: human and mouse at the moment")
  }
  #
  checkCategory(category)
  #
  if(!is_number(maxSize) | !is_number(minSize) | !is_number(counts) |!is_number(FDR)){
    stop("please provide numeric value")
  }
  #
  GOout <- list()
  #Extract all results
  #if(!is.null(ads)){
  #  results <- anota2seqGetDirectedRegulations(ads)
  #  #
  #  res <- vector("list", length = length(regulation))
  #  for(i in unique(contrast)){
  #    resTmp <- results[[i]][regulation[contrast==i]]
  #    res[which(contrast==i)] <- resTmp
  #  }
  #  names(res) <- paste(regulation, paste('c', contrast,sep=''), sep='_')
  #  if(!is.null(geneList)){
  #    res <- append(res, geneList)
  #  }
  #} else {
  #  res <- geneList
  #}
  #
  #if(!is.null(ads)){
  #  bg <- row.names(ads@dataP)
  #  if(!is.null(geneList)){
  #    bg <- unique(c(bg, as.character(unlist(geneList))))
  #  }
  #} else if(!is.null(customBg)){
  #  bg <- customBg
  #} else {
  #  stop("please provide background genes")
  #}
  res  <- a2sU_dataIn(a2sU)
  bg <- unlist(a2sU_bg(a2sU))
  if(length(setdiff(bg,unique(unlist(res))))==0){
    warning('Background seems not right as all genes are regulated')
  }
  
  #filter for slopes if indicated
  if(!is.null(genesSlopeFiltOut)){
    bg <- bg[!bg %in% genesSlopeFiltOut]
    res <- lapply(res, function(x) x[!x %in% genesSlopeFiltOut])
  }
  
  #Convert to enterezid
  bg_entrezID <- convertSymbolToEntrezID(geneList=bg, species=species)
  res_entrezID <- lapply(res,function(x) convertSymbolToEntrezID(geneList=x, species=species))
  
  GoLists <- res_entrezID[!sapply(res_entrezID, is.null)]

  #
  GOout  <- new("anota2seqUtilsGO",
              BP = NULL,
              CC = NULL,
              MF = NULL,
              KEGG = NULL)
  #
  for(sel in category){
    resOut <- list()
    if(sel=='KEGG'){
      for(i in 1:length(GoLists)){
        resTmp <- clusterProfiler::enrichKEGG(gene=GoLists[[i]],
                                        universe      = bg_entrezID,
                                        organism         = ifelse(species=="human","hsa",ifelse(species=="mouse","mmu",'ups')),
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 1,
                                        qvalueCutoff  = 1,
                                        minGSSize = minSize,
                                        maxGSSize = maxSize)
        #
        resOut[[names(GoLists)[i]]] <- resTmp
      }
    } else if(sel=='BP'|sel=='MF'|sel=='CC'){
      for(i in 1:length(GoLists)){
        resTmp <- clusterProfiler::enrichGO(gene=GoLists[[i]],
                                            universe      = bg_entrezID,
                                            OrgDb         = ifelse(species=="human","org.Hs.eg.db",ifelse(species=="mouse","org.Mm.eg.db",'ups')),
                                            ont           = sel,
                                            pAdjustMethod = "BH",
                                            pvalueCutoff  = 1,
                                            qvalueCutoff  = 1,
                                            minGSSize = minSize,
                                            maxGSSize = maxSize,
                                            readable = TRUE,
                                            pool = FALSE)
        resOut[[names(GoLists)[i]]] <- resTmp
      }
    }
    #remove counts below 10 and recalculate adjp and reformat
    for(i in 1:length(resOut)){
      print(i)
      tabTmp <- resOut[[i]]@result
      #
      tabTmp <- tabTmp[tabTmp$Count > counts,]
      tabTmp <- tabTmp[tabTmp$p.adjust < FDR, ]
      
      if(nrow(tabTmp)>0){
        #tabTmp$p.adjust <- stats::p.adjust( tabTmp$pvalue, method = 'BH')
        tabTmp$Size <- as.numeric(sub("\\/.*", "", tabTmp$BgRatio))
        tabTmp <- tabTmp[,c(1,2,9,10,5,6,8)]
      
        #tabTmp <- tabTmp[tabTmp$p.adjust < FDR, ]
        
        geneIDs_temp <- tabTmp$geneID
      
        checkID <- check_id_type(seqinr::c2s(strsplit(geneIDs_temp[[1]],'/')[[1]][1:5]))
        if(checkID=="entrezID"){
          tabTmp$geneID <- sapply(geneIDs_temp, function(x) paste(sort(convertEntrezIDToSymbol(unlist(strsplit(x,'/')),species=species)),collapse=':'),USE.NAMES = F)
        } else {
          tabTmp$geneID <- sapply(geneIDs_temp, function(x) paste(sort(unlist(strsplit(x,'/'))),collapse=':'),USE.NAMES = F)
        }
        #
        resOut[[i]]@result <- tabTmp
        resWrite <- lapply(resOut, function(x) x@result)
        #
        nameOut <- ifelse(is.null(name), paste('GO_',sel, '.xlsx',sep=''), paste(name,'_GO_', sel, '.xlsx',sep=''))
        WriteXLS::WriteXLS(resWrite,SheetNames = names(resWrite), ExcelFileName = nameOut, row.names=FALSE)
      } else {
        message('No significant results')
      }
    }
    slot(GOout, sel) <- resOut
  }
  a2sU@analysis@GO <- GOout
  #
  return(a2sU)
}

######
goDotplot <- function(a2sU,
                      #goIn,
                      category,
                      pool = TRUE,
                      termSel = NULL,
                      #colours = NULL,
                      nCategories=10,
                      size = 'Count',
                      pdfName=NULL){
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  if(!is_logical(pool)){
    stop("'pool' can only be only be logical: TRUE of FALSE ")
  }
  if(!is_number(nCategories)) {
    stop("please provide numeric value")
  }
  check_size(size)
  #
  checkCategory(category)
  #
  for(sel in category){
    if(is.null(slot(a2sU@analysis@GO, sel))){
      stop("Please run goAnalysis first for selected category")
    } else {
      goIn <- slot(a2sU@analysis@GO, sel)
    }
    #
    if(isTRUE(pool)){
      nameOut<- ifelse(is.null(pdfName), "pooled_GOdotplot.pdf", paste(pdfName, "pooled_GOdotplot.pdf", sep = "_"))
      #
      goDf <- data.table::rbindlist(lapply(goIn, function(x) x@result),use.names=TRUE, idcol=TRUE)
      colnames(goDf)[1] <- 'regulation'
      #
      if(!is.null(termSel)){
        goDf <- goDf[goDf$ID %in% termSel,]
      }
      #
      idx <- order(goDf$p.adjust, decreasing = FALSE)
      #
      goDf <- goDf[idx,]
      if(nCategories<nrow(goDf)){
        goDf <- goDf[1:nCategories,]
      }
    
      #to plot
      goDf$log10fdr <- -log10(goDf$p.adjust)
    
      #
      if(size=='geneRatio'){
        goDf$scale <- goDf$Count/goDf$Size
      } else {
        goDf$scale <- goDf$Count
      }

      #
      colOut <- colPlot(a2sU)[-1]
      names(colOut) <- names(a2sU@dataIn@geneList)
      
      pdf(nameOut, width = 8, height = 8, useDingbats = F)
      par(mar = c(5, 5, 3, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
      pOut <- ggplot2::ggplot(goDf, ggplot2::aes(x=log10fdr, y=reorder(Description,log10fdr), size=scale,colour = regulation)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(values = colOut) +   
        ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major = ggplot2::element_line(linetype = 'dashed', linewidth = 0.25), panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_blank(), legend.key.size = ggplot2::unit(0.5, 'cm')) +   
        ggplot2::xlab('-log10 FDR') +
        ggplot2::ylab(" ") 
        #+
        #ggplot2::ggtitle(paste(names(goIn)[i],sel, sep='_'))
      plot(pOut)
      dev.off()
    } else{
      for(i in 1:length(goIn)){
        nameOut <- ifelse(is.null(pdfName), paste(names(goIn)[i], "GOdotplot.pdf",sep='_'), paste(pdfName,names(goIn)[i],"GOdotplot.pdf", sep = "_"))
        goDf <- goIn[[i]]@result
        #
        if(!is.null(termSel)){
          goDf <- goDf[goDf$ID %in% termSel,]
        }
        #
        idx <- order(goDf$p.adjust, decreasing = FALSE)
        #
        goDf <- goDf[idx,]
        if(nCategories<nrow(goDf)){
          goDf <- goDf[1:nCategories,]
        }
        #to plot
        goDf$log10fdr <- -log10(goDf$p.adjust)
      
      #
        if(size=='geneRatio'){
          goDf$scale <- goDf$Count/goDf$Size
        } else {
          goDf$scale <- goDf$Count
        }
        
        #
        colOut <- colPlot(a2sU)[-1]
        names(colOut) <- names(a2sU_dataIn(a2sU))
        #
        pdf(nameOut, width = 8, height = 8, useDingbats = F)
        par(mar = c(5, 5, 3, 3), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.main = 1.7, cex.lab = 1)
        pOut <- ggplot2::ggplot(goDf, ggplot2::aes(x=log10fdr, y=reorder(Description,log10fdr), size=scale)) +
          ggplot2::geom_point(color= colOut[names(goIn)[i]]) +
          ggplot2::scale_color_manual(values = colOut[names(goIn)[i]]) +   
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid.major = ggplot2::element_line(linetype = 'dashed', linewidth = 0.25), panel.grid.minor = ggplot2::element_blank(),panel.background = ggplot2::element_blank(), legend.key.size = ggplot2::unit(0.5, 'cm')) +   
          ggplot2::xlab('-log10 FDR') +
          ggplot2::ylab(" ") +
          ggplot2::ggtitle(paste(names(goIn)[i],sel, sep='_')) 
        print(pOut)
        dev.off()
      }
    }
  }
}
