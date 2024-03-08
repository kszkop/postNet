gageAnalysis <- function(ads,
                         regulationGen,
                         contrastSel,
                         genesSlopeFiltOut=NULL,
                         rankIn = NULL,
                         species,
                         geneSetName
                         
){
  #
  gseaOut <- list()
  #
  if(!is.null(ads)){
    tmpAds <- anota2seq::anota2seqGetOutput(ads,
                                            analysis = regulationGen,
                                            output = "full",
                                            selContrast = contrastSel,
                                            getRVM = TRUE)
    #
    if (!is.null(genesSlopeFiltOut)) {
      tmpAdsFilt <- tmpAds[!row.names(tmpAds) %in% genesSlopeFiltOut, ]
    }  else {
      tmpAdsFilt <- tmpAds
    }
    #
    tmpP <- tmpAdsFilt[, "apvRvmP"]
    tmpEff <- tmpAdsFilt[, "apvEff"]
    #rankedRVMP <- rank(-log10(tmpP) * sign(tmpEff))
    rankIn <- tmpEff[order(tmpEff,decreasing = T)]
  } else if (!is.null(rankIn)){
    rankIn <- rankIn[order(rankIn,decreasing = T)]
  } else {
    stop("No anota2seq object or ranks provided")
  }
  #library(msigdb)
  #library(ExperimentHub)
  #library(AnnotationHub)
  #library(GSEABase)
  
  #library(fgsea)
  
  eh <- ExperimentHub::ExperimentHub()
  AnnotationHub::query(eh , 'msigdb')
  
  versionTmp <- as.character(sort(as.numeric(msigdb::getMsigdbVersions()),decreasing = T))[1]
  msigdbOut <- msigdb::getMsigdb(org = ifelse(species=="human",'hs', 'mm'), id = 'SYM', version = versionTmp)
  msigdbOut <- msigdb::appendKEGG(msigdbOut, version = versionTmp)
  #
  
  #Start with hallmark
  collectionTmp <- msigdb::subsetCollection(msigdbOut, collection = c('c2','h'))
  msigdb_ids = geneIds(collectionTmp)
  fgseaRes_hallmark <- fgsea(pathways = msigdb_hallmarks_ids, stat = fcOut, minSize  = 5, maxSize  = 500)
  head(fgseaRes_hallmark[order(padj), ])
  data.table::fwrite(fgseaRes_hallmark, file="smyd5_gsea_hallmark.txt", sep="\t", sep2=c("", " ", ""))
  
  
  c6_kegg <- subsetCollection(msigdb.hs, 'c6', 'CP:KEGG')
  msigdb_c6_kegg_ids = geneIds(c6_kegg)
  fgseaRes_c6_kegg<- fgsea(pathways = msigdb_c6_kegg_ids, stat = fcOut, minSize  = 5, maxSize  = 500)
  data.table::fwrite(fgseaRes_c6_kegg, file="smyd5_gsea_c6_kegg.txt", sep="\t", sep2=c("", " ", ""))
  
  c6 <- subsetCollection(msigdb.hs, 'c6')
  msigdb_c6_ids = geneIds(c6)
  fgseaRes_c6<- fgsea(pathways = msigdb_c6_ids, stat = fcOut, minSize  = 5, maxSize  = 500)
  
  outSet <- list()
  for(i in 1:length(selSet[-c(11,12)])){
    posOut <- which(grepl(selSet[-c(11,12)][i],names(msigdb.hs)))
    outSetTmp <- msigdb.hs[[posOut]]
    outSet[[i]] <- outSetTmp@geneIds
  }
  names(outSet) <- selSet[-c(11,12)]
  fgseaRes_selSet<- fgsea(pathways = outSet, stat = fcOut, minSize  = 5, maxSize  = 500)
  data.table::fwrite(fgseaRes_selSet, file="smyd5_gsea_selSet.txt", sep="\t", sep2=c("", " ", ""))
  
  ###Combine all sel
  toRunAll <- c(msigdb_hallmarks_ids,msigdb_c6_kegg_ids,outSet)
  fgseaRes_toRunAll <- fgsea(pathways = toRunAll, stat = fcOut, minSize  = 5, maxSize  = 500)
  data.table::fwrite(fgseaRes_toRunAll, file="smyd5_gsea_allCat.txt", sep="\t", sep2=c("", " ", ""))
  
  plotEnrichment <- function (pathway, stats, gseaParam = 1, ticksSize = 0.3) 
  {
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    g <- ggplot(toPlot, aes(x = x, y = y)) + geom_line(color = "firebrick1", linetype=1, size = 0.75)  + geom_line(color = "firebrick1") + theme_bw() + 
      theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + 
      theme(axis.line.x = element_line(color="black", linewidth = 0.5),
            axis.line.y = element_line(color="black", linewidth = 0.5)) +
      labs(x = "rank", y = "enrichment score") +
      geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
                                                                 y = -diff/4, xend = x, yend = diff/4), size = ticksSize, color="firebrick1")
    g
  }
  
  
  
  
  
}

