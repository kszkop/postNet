
setGeneric("a2sU_sequences",
           function(x, region) standardGeneric("a2sU_sequences"))
setMethod("a2sU_sequences", "anota2seqUtilsData",
          function(x, region){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5' or alternatively 'CCDS' if codon analysis performed with CCDS annotation.")
            }
            tmpReg <- slot(x@annot, region)
            seqOut <- tmpReg@seq
            return(seqOut)
          })

setGeneric("a2sU_id",
           function(x, region) standardGeneric("a2sU_id"))
setMethod("a2sU_id", "anota2seqUtilsData",
          function(x, region){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5' or alternatively 'CCDS' if codon analysis performed with CCDS annotation.")
            }
            tmpReg <- slot(x@annot, region)
            idOut <- tmpReg@id
            return(idOut)
          })

setGeneric("a2sU_geneID",
           function(x, region) standardGeneric("a2sU_geneID"))
setMethod("a2sU_geneID", "anota2seqUtilsData",
          function(x, region){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5' or alternatively 'CCDS' if codon analysis performed with CCDS annotation.")
            }
            tmpReg <- slot(x@annot, region)
            geneIDOut <- tmpReg@geneID
            return(geneIDOut)
          })

setGeneric("a2sU_dataIn",
           function(x) standardGeneric("a2sU_dataIn"))
setMethod("a2sU_dataIn", "anota2seqUtilsData",
          function(x){
            x@dataIn@geneList
          })

setGeneric("a2sU_bg",
           function(x) standardGeneric("a2sU_bg"))
setMethod("a2sU_bg", "anota2seqUtilsData",
          function(x){
            x@dataIn@background
          })

setGeneric("a2sU_eff",
           function(x) standardGeneric("a2sU_eff"))
setMethod("a2sU_eff", "anota2seqUtilsData",
          function(x){
            x@dataIn@effect
          })

setGeneric("a2sU_colours",
           function(x) standardGeneric("a2sU_colours"))
setMethod("a2sU_colours", "anota2seqUtilsData",
          function(x){
            x@dataIn@colours
          })

setGeneric("a2sU_species",
           function(x) standardGeneric("a2sU_species"))
setMethod("a2sU_species", "anota2seqUtilsData",
          function(x){
            x@species
          })

setGeneric("a2sU_version",
           function(x) standardGeneric("a2sU_version"))
setMethod("a2sU_version", "anota2seqUtilsData",
          function(x){
            x@version
          })

setGeneric("a2sU_selection",
           function(x) standardGeneric("a2sU_selection"))
setMethod("a2sU_selection", "anota2seqUtilsData",
          function(x){
            x@selection
          })

setGeneric("a2sU_motifs",
           function(x, region) standardGeneric("a2sU_motifs"))
setMethod("a2sU_motifs", "anota2seqUtilsData",
          function(x, region){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5'")
            }
            tmpReg <- slot(x@analysis@motifs, region)
            motifsOut <- tmpReg$motifsOut
            return(motifsOut)
          })

setGeneric("a2sU_motifsAnalysis",
           function(x, region, geneList) standardGeneric("a2sU_motifsAnalysis"))
setMethod("a2sU_motifsAnalysis", "anota2seqUtilsData",
          function(x, region, geneList){
            checkRegion(region)
            if(length(region)>1){
              stop("'region' can be only one of these: 'UTR3', 'CDS', 'UTR5'")
            }
            if(!any(geneList %in% names(a2sU_dataIn(x)))){
              stop('the regulatory geneList not in a2sU')
            }
            tmpReg <- slot(x@analysis@motifs, region)
            motifsAnalysisOut <- tmpReg[[geneList]]
            return(motifsAnalysisOut)
          })

setGeneric("a2sU_codonsAll",
           function(x) standardGeneric("a2sU_codonsAll"))
setMethod("a2sU_codonsAll", "anota2seqUtilsData",
          function(x){
            if(!checkUtils(x)){
              stop("It is not valid anota2seqUtils object")
            } else {
              tmpOut <- x@analysis@codons@codonsAll
              
              out <- s4_to_dataframe(tmpOut)
              return(out)
            }
          })

setGeneric("a2sU_codonsSel",
           function(x) standardGeneric("a2sU_codonsSel"))
setMethod("a2sU_codonsSel", "anota2seqUtilsData",
          function(x, comparison){
            if(!checkUtils(x)){
              stop("It is not valid anota2seqUtils object")
            } else {
              tmpOut <- x@analysis@codons@codons@
              
              out <- s4_to_dataframe(tmpOut)
              return(out)
            }
          })


a2sU_miRNA <- function(a2sU,
                       direction,
                       threshold) {
  #
  checkDirection(tolower(direction))
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  if(!is_number(threshold)){
    stop("'threshold' must be a number")
  }
  if(is.null(slot(a2sU@analysis, 'miRNA'))){
    stop("Please run mi first miRNAanalysis")
  } else {
    miRNAres <- a2sU@analysis@miRNA@miRNA_analysis
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


a2sU_miRNA_geneTargets <- function(a2sU,
                                   miRNAs){
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  if(is.null(slot(a2sU@analysis, 'miRNA_analysis'))){
    stop("Please run mi first miRNAanalysis")
  } else {
    miRNATmp <- a2sU@analysis@miRNA@miRNA_to_gene
  }
  miRNAsOut <- miRNATmp[which(names(miRNATmp) %in% miRNAs)]
  
  return(miRNAsOut)
}

###
a2sU_GO <- function(a2sU,
                    category,
                    geneList,
                    threshold) {
  #
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  checkCategory(category)
  if(!is_number(threshold)){
    stop("'threshold' must be a number")
  }
  #
  if(!any(geneList %in% names(a2sU_dataIn(a2sU)))){
    stop('None of the regulatory geneList in a2sU')
  }
  #
  if(is.null(slot(a2sU@analysis, 'GO'))){
    stop("Please run mi first GO analysis")
  } else {
    GOres <- slot(a2sU@analysis@GO,category)
  }
  GOresOut <- GOres[[which(geneList %in% names(GOres))]]@result
  GOresOut <- GOresOut[which(GOresOut[,6]< threshold),]

  if(nrow(GOresOut)>0){
    GOresOut <-  data.frame(GOresOut,row.names = NULL)
  } else {
    message('there are no GO categories to output')
  }
  return(GOresOut)
}

###
a2sU_gsea <- function(a2sU,
                      threshold=NULL) {
  
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  if(!is.null(threshold)){
    if(!is_number(threshold)){
      stop("'threshold' must be a number")
    }
  }
  if(is.null(slot(a2sU@analysis, 'GSEA'))){
    stop("Please run GSEA analysis first")
  } else {
    gseaOut <- slot(a2sU@analysis, 'GSEA')
  }
  if(!is.null(threshold)){
    gseaOut <- gseaOut[which(gseaOut[,8] < threshold),]
  }
  if(nrow(gseaOut)==0){
    message('there are no GSEA gene sets to output')
  }
  return(gseaOut)
}

a2sU_gage <- function(a2sU,
                      category,
                      direction,
                      threshold) {
  #
  checkDirection(tolower(direction))
  if(length(direction) != 1){
    stop("Please provide only one: greater or leas")
  }
  if (!checkUtils(a2sU)) {
    stop("a2sU is not a valid 'anota2seqUtilsData' object.")
  }
  if(!is_number(threshold)){
    stop("'threshold' must be a number")
  }
  checkCategory(category)
  if(length(category) != 1){
    stop("Please provide only one category")
  }
  if(is.null(slot(a2sU@analysis, 'GAGE'))){
    stop("Please run mi first miRNAanalysis")
  } else {
    GAGEres <- slot(a2sU@analysis@GAGE, category)
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

