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
      #for(pw in 1:length(pathwaysIn)){
      #  pathwaysIn[[pw]] <-  convertEntrezIDToSymbol(pathwaysIn[[pw]],species = species)
      #}
    }
    resOut <- gage::gage(exprs = rankIn,
                         gsets = pathwaysIn,
                         set.size= c(minSize, maxSize),
                         rank.test = TRUE,
                         use.fold = FALSE)
      
    slot(GAGEout, sel) <- resOut
  }
  a2sU@analysis@GAGE <- GAGEout
  #
  return(a2sU)
}

