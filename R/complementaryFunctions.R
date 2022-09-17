# This function assumes that anota2seqRegModes has been run on the ads already...
# This function is used to generate a list of genes for all contrasts and regulatory modes . . .
anota2seqGetDirectedRegulations <- function(ads){
  # Create a list for each contrast
  regModeList <- rep(list(NA),ncol(ads@contrasts))
  
  for(c in 1:ncol(ads@contrasts)){
    # Get the proper regulatory modes and store them in the list.
    translationUp <- ads@selectedTranslation@selectedRvmData[[c]][
      ads@selectedTranslation@selectedRvmData[[c]][,"apvEff"] > 0 &
        ads@selectedTranslation@selectedRvmData[[c]][,"singleRegMode"] == "translation",
      ]
    translationDown <-ads@selectedTranslation@selectedRvmData[[c]][
      ads@selectedTranslation@selectedRvmData[[c]][,"apvEff"] < 0 &
        ads@selectedTranslation@selectedRvmData[[c]][,"singleRegMode"] == "translation",
      ]
    translatedmRNAUp <- ads@selectedTranslatedmRNA@selectedRvmData[[c]][
      ads@selectedTranslatedmRNA@selectedRvmData[[c]][,"apvEff"] > 0 &
        ads@selectedTranslatedmRNA@selectedRvmData[[c]][,"singleRegMode"] == "translation",
      ]
    translatedmRNADown <-ads@selectedTranslatedmRNA@selectedRvmData[[c]][
      ads@selectedTranslatedmRNA@selectedRvmData[[c]][,"apvEff"] < 0 &
        ads@selectedTranslatedmRNA@selectedRvmData[[c]][,"singleRegMode"] == "translation",
      ]
    bufferingmRNAUp <- ads@selectedBuffering@selectedRvmData[[c]][
      ads@selectedBuffering@selectedRvmData[[c]][,"apvEff"] > 0 &
        ads@selectedBuffering@selectedRvmData[[c]][,"singleRegMode"] == "buffering",
      ]
    bufferingmRNADown <- ads@selectedBuffering@selectedRvmData[[c]][
      ads@selectedBuffering@selectedRvmData[[c]][,"apvEff"] < 0 &
        ads@selectedBuffering@selectedRvmData[[c]][,"singleRegMode"] == "buffering",
      ]
    mRNAAbundanceUp <- ads@mRNAAbundance@translatedmRNA[[c]][
      ads@mRNAAbundance@translatedmRNA[[c]][,"apvEff"] > 0 &
        ads@mRNAAbundance@translatedmRNA[[c]][,"singleRegMode"] == "abundance",
      ]
    mRNAAbundanceDown <- ads@mRNAAbundance@translatedmRNA[[c]][
      ads@mRNAAbundance@translatedmRNA[[c]][,"apvEff"] < 0 &
        ads@mRNAAbundance@translatedmRNA[[c]][,"singleRegMode"] == "abundance",
      ]
    totalmRNAUp <- ads@selectedTotalmRNA@selectedRvmData[[c]][
        ads@selectedTotalmRNA@selectedRvmData[[c]][,"apvEff"] > 0 &
          ads@selectedTotalmRNA@selectedRvmData[[c]][,"singleRegMode"] == "abundance",
      ]
    totalmRNADown <-ads@selectedTotalmRNA@selectedRvmData[[c]][
      ads@selectedTotalmRNA@selectedRvmData[[c]][,"apvEff"] < 0 &
        ads@selectedTotalmRNA@selectedRvmData[[c]][,"singleRegMode"] == "abundance",
      ]
      
      
    regModeList[[c]] <- list("translationUp" = rownames(translationUp),
                             "translationDown" = rownames(translationDown),
                             "translatedmRNAUp" = rownames(translatedmRNAUp),
                             "translatedmRNADown" = rownames(translatedmRNADown),
                             "bufferingmRNAUp" = rownames(bufferingmRNAUp),
                             "bufferingmRNADown" = rownames(bufferingmRNADown),
                             "mRNAAbundanceUp" = rownames(mRNAAbundanceUp),
                             "mRNAAbundanceDown" = rownames(mRNAAbundanceDown),
                             "totalmRNAUp" = rownames(totalmRNAUp),
                             "totalmRNADown" = rownames(totalmRNADown))#,
                             #"background" = setdiff(rownames(ads@dataP),c(rownames(translationUp),
                                                                          #rownames(translationDown),
                                                                          #rownames(bufferingmRNADown),
                                                                          #rownames(bufferingmRNAUp),
                                                                          #rownames(mRNAAbundanceUp),
                                                                          #rownames(mRNAAbundanceDown))
                             #))
    
  }
  
  return(regModeList)
  
}


###
checkAvailableVersions <- function(species){
  #list existing species
  currTmp <- list.files(system.file("extdata/annotation/refseq",package = "anota2seqUtils"))
  
  if(!species %in% currTmp){
    stop("The pre-prepared annotation file for that species does not exist. Please use option createFromFile")
  }  else {
    dir <- system.file(paste("extdata/annotation/refseq", species,sep='/'),package = "anota2seqUtils")
    listVersions <- list.files(path = dir)
    print(listVersions)
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
calc_motif <- function(x, motif, len){
  seqTmp <- x
  motOut  <- seqinr::words.pos(motif,seqTmp)
  if(length(motOut)>0){
    #Check overlapping and collapse them
    gROut <- GenomicRanges::reduce(GenomicRanges::GRanges(seqnames='tmp', ranges=IRanges::IRanges(start=motOut,end=motOut)),min.gapwidth=len)
    #
    nMot <- length(gROut)
  } else {
    nMot <- 0
  }
  return(nMot)
}

#Combine 
calc_g4 <- function(x,min_score){
  seqTmp <- DNAString(x)
  predTmp <- pqsfinder::pqsfinder(seqTmp, min_score = min_score, strand = '+')
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
codonCount <- function(gene,seq){
  #
  tmpEff <- seqinr::uco(seqinr::s2c(seq),index = "eff")
  tmpFreq <- seqinr::uco(seqinr::s2c(seq),index = "freq")
  
  tmpCodon <- data.frame(geneID=gene,codon=toupper(names(tmpEff)),AA=seqinr::aaa(seqinr::translate(seqinr::s2c(seqinr::c2s(toupper(names(tmpEff)))))), codonCount=as.numeric(tmpEff),codonFreq=as.numeric(tmpFreq))
  
  tmpCodon <- tmpCodon %>% group_by(AA) %>% mutate(AACountPerGene=sum(codonCount))
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
  
    for(AA in 1:length(uniqAA)){
      #
      tmpDf <- df[df$AA == uniqAA[AA],]
      #
      codons <- as.character(df$codon[df$AA == uniqAA[AA]])
    
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

combSeq <- function(seq1, seq2){
  seq1_tmp <- seqinr::s2c(seq1)
  seq2_tmp <- seqinr::s2c(seq2)
  seqOut <- seqinr::c2s(c(seq1_tmp,seq2_tmp))
  return(seqOut)
}

calc_uORF <- function(seqTmp, ext, context){
  #
  nTmp <- as.numeric()
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
        stopOut <- stopOut[which((stopOut - startOut[i])>0)]
        #
        potORF <- stopOut - startOut[i]
        #check in frame 
        inFrameCheck <- potORF %% 3
        #Take first stop in frame
        stopOut <- stopOut[which(inFrameCheck==0)]
        #if maybe postition needed
        #if(length(stopOut)>0){
        #  stopOut <- min(stopOut)+2
        #}
        #
        if(length(stopOut)>0){
          nTmp[i] <- startOut[i] 
        }
      }
      #
      nOut <- length(nTmp)
    } else {
      #
      nOut <- 0
    }
  } else {
    nOut <- 0
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
