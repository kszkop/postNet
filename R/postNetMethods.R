setGeneric("ptn_sequences",
           function(x, region) standardGeneric("ptn_sequences"))
setMethod("ptn_sequences", "postNetData",
          function(x, region){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5' or alternatively 'CCDS' if codon analysis performed with CCDS annotation.")
            }
            tmpReg <- slot(x@annot, region)
            seqOut <- tmpReg@sequences
            return(seqOut)
          })

setGeneric("ptn_id",
           function(x, region) standardGeneric("ptn_id"))
setMethod("ptn_id", "postNetData",
          function(x, region){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5' or alternatively 'CCDS' if codon analysis performed with CCDS annotation.")
            }
            tmpReg <- slot(x@annot, region)
            idOut <- tmpReg@id
            return(idOut)
          })

setGeneric("ptn_geneID",
           function(x, region) standardGeneric("ptn_geneID"))
setMethod("ptn_geneID", "postNetData",
          function(x, region){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5' or alternatively 'CCDS' if codon analysis performed with CCDS annotation.")
            }
            tmpReg <- slot(x@annot, region)
            geneIDOut <- tmpReg@geneID
            return(geneIDOut)
          })

setGeneric("ptn_dataIn",
           function(x) standardGeneric("ptn_dataIn"))
setMethod("ptn_dataIn", "postNetData",
          function(x){
            x@dataIn
          })

setGeneric("ptn_geneList",
           function(x) standardGeneric("ptn_geneList"))
setMethod("ptn_geneList", "postNetData",
          function(x){
            x@dataIn@geneList
          })

setGeneric("ptn_background",
           function(x) standardGeneric("ptn_background"))
setMethod("ptn_background", "postNetData",
          function(x){
            x@dataIn@background
          })

setGeneric("ptn_effect",
           function(x) standardGeneric("ptn_effect"))
setMethod("ptn_effect", "postNetData",
          function(x){
            x@dataIn@effect
          })

setGeneric("ptn_colours",
           function(x) standardGeneric("ptn_colours"))
setMethod("ptn_colours", "postNetData",
          function(x){
            x@dataIn@colours
          })

setGeneric("ptn_species",
           function(x) standardGeneric("ptn_species"))
setMethod("ptn_species", "postNetData",
          function(x){
            x@species
          })

setGeneric("ptn_version",
           function(x) standardGeneric("ptn_version"))
setMethod("ptn_version", "postNetData",
          function(x){
            x@version
          })

setGeneric("ptn_selection",
           function(x) standardGeneric("ptn_selection"))
setMethod("ptn_selection", "postNetData",
          function(x){
            x@selection
          })

setGeneric("ptn_motifSelection",
           function(x, region) standardGeneric("ptn_motifSelection"))
setMethod("ptn_motifSelection", "postNetData",
          function(x, region){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5'")
            }
            tmpReg <- slot(x@analysis@motifs, region)
            motifsOut <- tmpReg$motifSelection
            return(motifsOut)
          })

setGeneric("ptn_motifgeneList",
           function(x, region, geneList) standardGeneric("ptn_motifgeneList"))
setMethod("ptn_motifgeneList", "postNetData",
          function(x, region, geneList){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5'")
            }
            if(!any(geneList %in% names(ptn_dataIn(x)))){
              stop('the regulatory geneList not in ptn')
            }
            tmpReg <- slot(x@analysis@motifs, region)
            motifsAnalysisOut <- tmpReg[[geneList]]
            return(motifsAnalysisOut)
          })

setGeneric("ptn_codonAnalysis",
           function(x) standardGeneric("ptn_codonAnalysis"))
setMethod("ptn_codonAnalysis", "postNetData",
          function(x){
            if(!checkPtn(x)){
              stop("It is not valid postNet object")
            } else {
              tmpOut <- x@analysis@codons@codonAnalysis
              
              out <- s4_to_dataframe(tmpOut)
              return(out)
            }
          })

setGeneric("ptn_codonSelection",
           function(x,comparison) standardGeneric("ptn_codonSelection"))
setMethod("ptn_codonSelection", "postNetData",
          function(x, comparison){
            if(!checkPtn(x)){
              stop("It is not valid postNet object")
            } else {
              out <- x@analysis@codons@codonSelection[[comparison]]

              return(out)
            }
          })

setGeneric("ptn_features",
           function(x) standardGeneric("ptn_features"))
setMethod("ptn_features", "postNetData",
          function(x){
            if(!checkPtn(x)){
              stop("It is not valid postNet object")
            } else {
              out <- x@features
              return(out)
            }
          })


ptn_miRNA_analysis <- function(ptn,
                               direction,
                               threshold) {
  #
  checkDirection(tolower(direction))
  if (!checkPtn(ptn)) {
    stop("ptn is not a valid 'postNetData' object.")
  }
  if(!is_number(threshold)){
    stop("'threshold' must be a number")
  }
  if(is.null(slot(ptn@analysis, 'miRNA'))){
    stop("Please run mi first miRNAanalysis")
  } else {
    miRNAres <- ptn@analysis@miRNA@miRNA_analysis
  }
  
  if(tolower(direction)=='greater'){
    resOut <- miRNAres$greater
  } else if (tolower(direction)=='less') {
    resOut <- miRNAres$less
  }
  #
  resOut <- resOut[which(resOut[,4]< threshold),]
  #
  if(nrow(resOut)>0){
    resOut <- resOut[,c(1,2,5,3,4)]
    resOut <-  data.frame(id=row.names(resOut),resOut,row.names = NULL)
  } else {
    message('there are no output miRNAs')
  }
  return(resOut)
}


ptn_miRNA_to_gene <- function(ptn,
                              miRNAs){
  #
  if (!checkPtn(ptn)) {
    stop("ptn is not a valid 'postNetData' object.")
  }
  if(is.null(slot(ptn@analysis, 'miRNA_analysis'))){
    stop("Please run mi first miRNAanalysis")
  } else {
    miRNATmp <- ptn@analysis@miRNA@miRNA_to_gene
  }
  miRNAsOut <- miRNATmp[which(names(miRNATmp) %in% miRNAs)]
  
  return(miRNAsOut)
}

###
ptn_GO <- function(ptn,
                    category,
                    geneList,
                    threshold) {
  #
  if (!checkPtn(ptn)) {
    stop("ptn is not a valid 'postNetData' object.")
  }
  checkCategory(category)
  if(!is_number(threshold)){
    stop("'threshold' must be a number")
  }
  #
  if(!any(geneList %in% names(ptn_dataIn(ptn)))){
    stop('None of the regulatory geneList in ptn')
  }
  #
  if(is.null(slot(ptn@analysis, 'GO'))){
    stop("Please run mi first GO analysis")
  } else {
    GOres <- slot(ptn@analysis@GO,category)
  }
  GOresOut <- GOres[[which(geneList == names(GOres))]]@result
  GOresOut <- GOresOut[which(GOresOut[,6]< threshold),]

  if(nrow(GOresOut)>0){
    GOresOut <-  data.frame(GOresOut,row.names = NULL)
  } else {
    message('there are no GO categories to output')
  }
  return(GOresOut)
}

###
ptn_GSEA<- function(ptn,
                      threshold=NULL) {
  
  if (!checkPtn(ptn)) {
    stop("ptn is not a valid 'postNetData' object.")
  }
  if(!is.null(threshold)){
    if(!is_number(threshold)){
      stop("'threshold' must be a number")
    }
  }
  if(is.null(slot(ptn@analysis, 'GSEA'))){
    stop("Please run GSEA analysis first")
  } else {
    gseaOut <- slot(ptn@analysis, 'GSEA')
  }
  if(!is.null(threshold)){
    gseaOut <- gseaOut[which(gseaOut[,8] < threshold),]
  }
  if(nrow(gseaOut)==0){
    message('there are no GSEA gene sets to output')
  }
  return(gseaOut)
}

ptn_GAGE <- function(ptn,
                      category,
                      direction,
                      threshold) {
  #
  checkDirection(tolower(direction))
  if(length(direction) != 1){
    stop("Please provide only one: greater or leas")
  }
  if (!checkPtn(ptn)) {
    stop("ptn is not a valid 'postNetData' object.")
  }
  if(!is_number(threshold)){
    stop("'threshold' must be a number")
  }
  checkCategory(category)
  if(length(category) != 1){
    stop("Please provide only one category")
  }
  if(is.null(slot(ptn@analysis, 'GAGE'))){
    stop("Please run mi first miRNAanalysis")
  } else {
    GAGEres <- slot(ptn@analysis@GAGE, category)
  }
  
  if(tolower(direction)=='greater'){
    resOut <- GAGEres$greater
  } else if (tolower(direction)=='less') {
    resOut <- GAGEres$less
  }
  #
  resOut <- resOut[which(resOut[,6] < threshold),]
  #
  if(nrow(resOut)==0){
    message('there are no enriched terms')
  } else {
    return(resOut)
  }
}

ptn_model <- function(ptn, analysis_type, model, comparison){
  if (!checkPtn(ptn)) {
    stop("ptn is not a valid 'postNetData' object.")
  }
  if(!is_valid_analysis_type(analysis_type)){
    stop("'analysis_type' can be only 'lm' for linear model or 'rf' for random forest")
  }
  if(!check_model(model, analysis_type = analysis_type)){
    stop("please provide correct model for a analysis type")
  }
  if (!is_number(comparison)) {
    stop("please provide correct comparison number")
    if(length(comparison) != 1){
      stop("comparison can be only one")
    }
  }
  tmpIn <- slot(ptn@analysis@featureIntegration,analysis_type)
  tmpIn <- tmpIn[[comparison]]
  #
  tmpOut <- slot(tmpIn,model)
  #
  return(tmpOut)
}


ptn_selectedFeatures <- function(ptn, analysis_type){
  if (!checkPtn(ptn)) {
    stop("ptn is not a valid 'postNetData' object.")
  }
  if(!is_valid_analysis_type(analysis_type)){
    stop("'analysis_type' can be only 'lm' for linear model or 'rf' for random forest")
  }
  
  tmpIn <- slot(ptn@analysis@featureIntegration,analysis_type)
  #
  tmpOut <- slot(tmpIn, "selectedFeatures")
  #
  return(tmpOut)
}

setGeneric("ptn_networkGraph",
           function(x) standardGeneric("ptn_networkGraph"))
setMethod("ptn_networkGraph", "postNetData",
          function(x){
            if(!checkPtn(x)){
              stop("It is not valid postNet object")
            } else {
              tmpOut <- x@analysis@featureIntegration@lm@networkGraph
              return(tmpOut)
            }
          })
