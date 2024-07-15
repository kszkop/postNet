gageAnalysis <- function(a2sU,
                         category, #To choose from (gene ontologies) 'BP', 'CC', 'MF' or 'KEGG'
                         genesSlopeFiltOut=NULL,
                         maxSize = 500,
                         minSize = 10
){
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  species <- a2sU_species(a2sU)
  if (!species %in% c("human","mouse")) {
    stop("This option is only  available for species: human and mouse at the moment")
  }
  if(!is_number(maxSize) | !is_number(minSize)){
    stop("please provide numeric value")
  }
  #
  gageOut <- list()

  effTmp <- a2sU_eff(a2sU)
  if (!is.null(genesSlopeFiltOut)) {
    effIn <- effTmp[!names(effTmp) %in% genesSlopeFiltOut ]
  }  else {
    effIn <- effTmp
  }
  #
  convTab <- convertSymbolToEntrezID(names(effIn),species=species)
  names(convTab) <- convertEntrezIDToSymbol(convTab,species=species)
  #
  names(effIn) <- as.character(convTab)[match(names(effIn), names(convTab))]
  #remove NA
  effIn <- effIn[!is.na(names(effIn))]
  
  #
  rankIn <- effIn[order(effIn,decreasing = T)]
  
  #
  GAGEout  <- new("anota2seqUtilsGAGE",
                BP = NULL,
                CC = NULL,
                MF = NULL,
                KEGG = NULL)
  
  
  for(sel in category){
    resOut <- list()
    if(sel=='KEGG'){
      pathwayKegg<- gage::kegg.gsets(species=species,id.type = "entrez")
      pathwaysIn <- pathwayKegg$kg.sets
      #for(pw in 1:length(pathwaysIn)){
      #  pathwaysIn[[pw]] <- convertEntrezIDToSymbol(pathwaysIn[[pw]],species = species)
      #}
    } else {
      if(!exists('pathwayGO')){
        pathwayGO <- gage::go.gsets(species=species,id.type = "EG")
      }
      if(sel == 'BP') {
        pathwaysIn <- pathwayGO$go.sets[pathwayGO$go.subs$BP]
      } else if(sel == 'MF') {
        pathwaysIn <- pathwayGO$go.sets[pathwayGO$go.subs$MF]
      } else if(sel == 'CC') {
        pathwaysIn <- pathwayGO$go.sets[pathwayGO$go.subs$CC]
      } 
    }
    resOut <- gage::gage(exprs = rankIn,
                         gsets = pathwaysIn,
                         set.size= c(minSize, maxSize),
                         rank.test = TRUE,
                         use.fold = FALSE)
    
    #
    pathwaysGenes <- list()
    for(pw in 1:length(pathwaysIn)){
      glist <- unlist(pathwaysIn[pw])
      tmpGenes <- glist[glist %in% names(rankIn)]
      if(length(tmpGenes)>0){
        tmpGenes_conv <- names(convTab)[match(as.character(tmpGenes),convTab)]
        #tmpGenes_conv <- convertEntrezIDToSymbol(pathwaysIn[[pw]], species = species)
        #AnnotationDbi::mapIds(x = org.Hs.eg.db, keys = as.character(glist), column="SYMBOL", keytype="ENTREZID")
      } else {
        tmpGenes_conv <- NA
      }
      pathwaysGenes[[pw]] <-  paste(tmpGenes_conv,collapse=':')
    }
    names(pathwaysGenes) <- names(pathwaysIn)
    pathwaysGenes <- pathwaysGenes[as.numeric(which(unlist(lapply(pathwaysGenes, function(x) x!='NA'))))]

    if(nrow(resOut$greater)>0){
      grTmpG <- resOut$greater
      grTmpG <- grTmpG[,c(1,2,5,3,4)]
      grTmpG <-  data.frame(id=row.names(grTmpG),grTmpG,row.names = NULL)
      
      grTmpG$Genes <- as.character(pathwaysGenes)[match(grTmpG$id, names(pathwaysGenes))]
      resOut$greater <- grTmpG
    }
    if(nrow(resOut$less)>0){
      grTmpL <- resOut$less
      grTmpL <- grTmpL[,c(1,2,5,3,4)]
      grTmpL <-  data.frame(id=row.names(grTmpL),grTmpL,row.names = NULL)
      
      grTmpL$Genes <- as.character(pathwaysGenes)[match(grTmpL$id, names(pathwaysGenes))]
      resOut$less <- grTmpL
    }
    #
    slot(GAGEout, sel) <- resOut
  }
  a2sU@analysis@GAGE <- GAGEout
  #
  return(a2sU)
}

