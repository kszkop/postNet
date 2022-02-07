goAnalysis <- function(ads,
                       regulation,
                       contrast,
                       genesSlopeFiltOut=NULL,
                       species,
                       category, #To choose from (gene ontologies) 'BP', 'CC', 'MF'  or 'KEGG'
                       maxSize=500,
                       minSize=10,
                       counts=10,
                       FDR=0.15,
                       myCond=F,
                       name
){
  #Subset annot for only expressed genes
  bg <- row.names(ads@dataP)
  
  #filter for slopes if indicated
  if(!is.null(genesSlopeFiltOut)){
    bg <- bg[!bg %in% genesSlopeFiltOut]
  }
  
  #Convert to enterezid
  bg_entrezID <- convertSymbolToEntrezID(geneList=bg, species=species)
  
  #Extract all results
  results <- anota2seqGetDirectedRegulations(ads)
  
  #
  res <- vector("list", length = length(regulation))
  for(i in 1:length(unique(contrast))){
    resTmp <- results[[i]][names(results[[i]]) %in% regulation[which(contrast==i)]]
    res[which(contrast==i)] <- resTmp
  }
  
  #There should not be any but just in case filter for slopes if indicated
  if(!is.null(genesSlopeFiltOut)){
    res <- lapply(res, function(x) x[!x %in% genesSlopeFiltOut])
  }
  
  #convert to entrezID
  res_entrezID <- list()
  for(i in 1:length(res)){
    symTmp <- res[[i]]
    res_entrezID[[i]] <- convertSymbolToEntrezID(geneList=symTmp, species=species)
  }
  names(res_entrezID) <- regulation
                       
  #genes <- res
  entrezIDList <- list()
  nullVec <- vector("numeric")                    
  
  for( i in 1:length(res_entrezID)){
    if(!is.null(res_entrezID[[i]])){
      nullVec <- c(nullVec,i)
    }
    entrezIDList[[i]] <- res_entrezID[[i]]
  }
  
  #genes <- genes[nullVec]
  #print(length(genes))
  names(entrezIDList) <- regulation
  #print(length(entrezIDList))
  geneBackground <- bg_entrezID
  
  GoLists <- entrezIDList[!sapply(entrezIDList, is.null)]
  
  ####
  if(category=='KEGG'){
    GOobjects <- list()
    for(i in 1:length(GoLists)){
      GOobjects[[i]]<- new("KEGGHyperGParams",
                           geneIds=as.integer(GoLists[[i]]), ### Genes from the 
                           universeGeneIds=as.integer(geneBackground),
                           annotation=ifelse(species=="human","org.Hs.eg.db",ifelse(species=="mouse","org.Mm.eg.db",'ups')),
                           #ontology=category,
                           pvalueCutoff=1,
                           #conditional=myCond,
                           testDirection="over") 
    }
  } else if(category=='BP'|category=='MF'|category=='CC'){
    #
    GOobjects <- list()
    for(i in 1:length(GoLists)){
      GOobjects[[i]]<- new("GOHyperGParams",
                         geneIds=as.integer(GoLists[[i]]), ### Genes from the 
                         universeGeneIds=as.integer(geneBackground),
                         annotation=ifelse(species=="human","org.Hs.eg.db",ifelse(species=="mouse","org.Mm.eg.db",'ups')),
                         ontology=category,
                         pvalueCutoff=1,
                         conditional=myCond,
                         testDirection="over") 
    }
  } else {
    stop('Wrong category! Please select from: BP, MF, CC, KEGG or PFAM')
  }
  #
  hgOver <- list()
  for(i in 1:length(GOobjects)){
    hgOver[[i]] <- GOstats::hyperGTest(GOobjects[[i]])
  }
  
  names(hgOver) <- names(GoLists)
  
  GOtables <- list()
  for( i in 1:length(hgOver)){
    GOtables[[i]] <- summary(hgOver[[i]])
  }
  names(GOtables)<- names(hgOver)
  print('KS1')
  #for( i in 1:length(GOtables)){
  #  rownames(GOtables[[i]])<- GOtables[[i]][,"Term"]
  #}
  
  ### Filter go tables output 
  ### Size of the GO terms
  ### number of genes per GO term
  ### OddsRatio
  ###

  for(i in 1:length(GOtables)){
    GOtables[[i]] <- GOtables[[i]][which(GOtables[[i]]["Size"] <= maxSize & GOtables[[i]]["Size"] > minSize& GOtables[[i]]["Count"] > counts),]
  }

  
  print('KS2')
  ### Calculate FDRs for the output
  for ( i in 1:length(hgOver)){
    
    FDR<- p.adjust(GOtables[[i]][,"Pvalue"],method="BH")
    GOtables[[i]]<- cbind(GOtables[[i]],FDR)
  }
  
  print('KS3')
  # ### Add genes to GO terms...
  genesPerTerm <- list()
  
  for(i in 1:length(hgOver)){
    genesPerTerm[[i]] <- Category::geneIdsByCategory(hgOver[[i]])
  }
  for( i in 1:length(genesPerTerm)){
    if(category=='KEGG'){
      genesPerTerm[[i]] <- genesPerTerm[[i]][GOtables[[i]][,"KEGGID"]]
    } else {
      genesPerTerm[[i]] <- genesPerTerm[[i]][GOtables[[i]][,paste("GO",category,"ID",sep='')]]
    }
  }
  print('KS4')
  #
  if(species=="human"){
    symbol2id <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL2EG)
  } else if(species=="mouse"){
    symbol2id <- as.list(org.Mm.eg.db::org.Mm.egSYMBOL2EG)
  }
  #
  id2symb <- as.list(setNames(names(symbol2id), symbol2id))
  print('KS5')
  geneNamesPerTerm <- genesPerTerm
  for(k  in 1:length(genesPerTerm)){
    for(l in 1:length(genesPerTerm[[k]])){
      if(length(genesPerTerm[[k]]) >0){
        tmpList <- as.vector(genesPerTerm[[1]][[1]],"character")
        geneNamesPerTerm[[k]][[l]]<- paste(sort(convertEntrezIDToSymbol(tmpList, species=species)),collapse=":")
      }
    }
  }
  print('KS6')
  #
  for(i in 1:length(GOtables)){
    genes <- unlist(geneNamesPerTerm[[i]])
    tmpTab <- cbind(GOtables[[i]],genes)
    #reformat
    tmpTab <- tmpTab[,c(1,7,2,8,3:6,9)]
    #row.names(tmpTab) <- NULL
    GOtables[[i]] <- tmpTab 
  }
  print('KS7')
  WriteXLS::WriteXLS(GOtables,SheetNames = names(GOtables),ExcelFileName = paste(name, paste(category, '.xlsx',sep=''),sep='_'),row.names=FALSE)
  return(list(GOtables,hgOver))
}


