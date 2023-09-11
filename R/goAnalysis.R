goAnalysis <- function(ads,
                       regulation,
                       contrast,
                       genesSlopeFiltOut=NULL,
                       geneList = NULL,
                       customBg = NULL,
                       species,
                       category, #To choose from (gene ontologies) 'BP', 'CC', 'MF'  or 'KEGG'
                       maxSize=500,
                       minSize=10,
                       counts=10,
                       FDR=0.15,
                       myCond=F,
                       name
){
  GOout <- list()
  #Extract all results
  if(!is.null(ads)){
    results <- anota2seqGetDirectedRegulations(ads)
    #
    res <- vector("list", length = length(regulation))
    for(i in unique(contrast)){
      resTmp <- results[[i]][regulation[contrast==i]]
      res[which(contrast==i)] <- resTmp
    }
    names(res) <- paste(regulation, paste('c', contrast,sep=''), sep='_')
    if(!is.null(geneList)){
      res <- append(res, geneList)
    }
  } else {
    res <- geneList
  }

  #
  if(!is.null(ads)){
    bg <- row.names(ads@dataP)
    if(!is.null(geneList)){
      bg <- unique(c(bg, as.character(unlist(geneList))))
    }
  } else if(!is.null(customBg)){
    bg <- customBg
  } else {
    stop("please provide background genes")
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
  
  ####
  for(j in 1:length(category)){
    categoryIn <- category[j]
    #
    GOobjects <- list()
    if(categoryIn=='KEGG'){
        for(i in 1:length(GoLists)){
          GOobjects[[i]]<- new("KEGGHyperGParams",
                           geneIds=as.integer(GoLists[[i]]), ### Genes from the 
                           universeGeneIds=as.integer(bg_entrezID),
                           annotation=ifelse(species=="human","org.Hs.eg.db",ifelse(species=="mouse","org.Mm.eg.db",'ups')),
                           pvalueCutoff=1,
                           testDirection="over") 
        }
    } else if(categoryIn=='BP'|categoryIn=='MF'|categoryIn=='CC'){
      for(i in 1:length(GoLists)){
        GOobjects[[i]]<- new("GOHyperGParams",
                             geneIds=as.integer(GoLists[[i]]), ### Genes from the 
                             universeGeneIds=as.integer(bg_entrezID),
                             annotation=ifelse(species=="human","org.Hs.eg.db",ifelse(species=="mouse","org.Mm.eg.db",'ups')),
                             ontology=categoryIn,
                             pvalueCutoff=1,
                             conditional=myCond,
                             testDirection="over") 
      }
    } else {
      stop('Wrong category! Please select from: BP, MF, CC, KEGG')
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
        
    #
    for(i in 1:length(GOtables)){
        GOtables[[i]] <- GOtables[[i]][which(GOtables[[i]][,"Size"] <= maxSize & GOtables[[i]][,"Size"] > minSize & GOtables[[i]][,"Count"] > counts),]
    }
    ### Calculate FDRs for the output
    for ( i in 1:length(hgOver)){
        FDR<- p.adjust(GOtables[[i]][,"Pvalue"],method="BH")
        GOtables[[i]]<- cbind(GOtables[[i]],FDR)
    }
        
    # ### Add genes to GO terms...
    genesPerTerm <- list()
    for(i in 1:length(hgOver)){
        genesPerTerm[[i]] <- Category::geneIdsByCategory(hgOver[[i]])
    }
    for( i in 1:length(genesPerTerm)){
      if(categoryIn=='KEGG'){
        genesPerTerm[[i]] <- genesPerTerm[[i]][GOtables[[i]][,"KEGGID"]]
      } else {
        genesPerTerm[[i]] <- genesPerTerm[[i]][GOtables[[i]][,paste("GO",categoryIn,"ID",sep='')]]
      }
    }
    #
    if(species=="human"){
        symbol2id <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL2EG)
    } else if(species=="mouse"){
        symbol2id <- as.list(org.Mm.eg.db::org.Mm.egSYMBOL2EG)
    }
    #
    id2symb <- as.list(setNames(names(symbol2id), symbol2id))

    geneNamesPerTerm <- genesPerTerm
    for(k  in 1:length(genesPerTerm)){
        for(l in 1:length(genesPerTerm[[k]])){
            if(length(genesPerTerm[[k]]) >0){
                tmpList <- as.vector(genesPerTerm[[1]][[1]],"character")
                geneNamesPerTerm[[k]][[l]]<- paste(sort(convertEntrezIDToSymbol(tmpList, species=species)),collapse=":")
            }
        }
    }
    #
    for(i in 1:length(GOtables)){
        genes <- unlist(geneNamesPerTerm[[i]])
        tmpTab <- cbind(GOtables[[i]],genes)
        #reformat
        tmpTab <- tmpTab[,c(1,7,2,8,3:6,9)]
        #row.names(tmpTab) <- NULL
        GOtables[[i]] <- tmpTab 
    }
    WriteXLS::WriteXLS(GOtables,SheetNames = names(GOtables),ExcelFileName = paste(name, paste(categoryIn, '.xlsx',sep=''),sep='_'),row.names=FALSE)
    GOcat <- list(GOtables,hgOver)
    GOout[[categoryIn]] <- GOcat
  }
  return(GOout)
}


