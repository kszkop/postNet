
# accessor regions data
setGeneric("anota2seqUtilsGetUTR5Seq",
           function(x) standardGeneric("anota2seqUtilsGetUTR5Seq"))
setMethod("anota2seqUtilsGetUTR5Seq", "anota2seqUtilsData",
            function(x){
              x@annot@UTR5@seq
            })

setGeneric("anota2seqUtilsGetUTR5id",
           function(x) standardGeneric("anota2seqUtilsGetUTR5id"))
setMethod("anota2seqUtilsGetUTR5id", "anota2seqUtilsData",
          function(x){
            x@annot@UTR5@id
          })

setGeneric("anota2seqUtilsGetUTR5geneID",
           function(x) standardGeneric("anota2seqUtilsGetUTR5geneID"))
setMethod("anota2seqUtilsGetUTR5geneID", "anota2seqUtilsData",
          function(x){
            x@annot@UTR5@geneID
          })

setGeneric("anota2seqUtilsGetCDSSeq",
           function(x) standardGeneric("anota2seqUtilsGetCDSSeq"))
setMethod("anota2seqUtilsGetCDSSeq", "anota2seqUtilsData",
          function(x){
            x@annot@CDS@seq
          })

setGeneric("anota2seqUtilsGetCDSid",
           function(x) standardGeneric("anota2seqUtilsGetCDSid"))
setMethod("anota2seqUtilsGetCDSid", "anota2seqUtilsData",
          function(x){
            x@annot@CDS@id
          })

setGeneric("anota2seqUtilsGetCDSgeneID",
           function(x) standardGeneric("anota2seqUtilsGetCDSgeneID"))
setMethod("anota2seqUtilsGetCDSgeneID", "anota2seqUtilsData",
          function(x){
            x@annot@CDS@geneID
          })

setGeneric("anota2seqUtilsGetUTR3Seq",
           function(x) standardGeneric("anota2seqUtilsGetUTR3Seq"))
setMethod("anota2seqUtilsGetUTR3Seq", "anota2seqUtilsData",
          function(x){
            x@annot@UTR3@seq
          })

setGeneric("anota2seqUtilsGetUTR3id",
           function(x) standardGeneric("anota2seqUtilsGetUTR3id"))
setMethod("anota2seqUtilsGetUTR3id", "anota2seqUtilsData",
          function(x){
            x@annot@UTR3@id
          })

setGeneric("anota2seqUtilsGetUTR3geneID",
           function(x) standardGeneric("anota2seqUtilsGetUTR3geneID"))
setMethod("anota2seqUtilsGetUTR3geneID", "anota2seqUtilsData",
          function(x){
            x@annot@UTR3@geneID
          })

setGeneric("anota2seqUtilsGetUTR3",
           function(x) standardGeneric("anota2seqUtilsGetUTR3"))
setMethod("anota2seqUtilsGetUTR3", "anota2seqUtilsData",
          function(x){
            x@annot@UTR3
          })

setGeneric("anota2seqUtilsGetUTR5",
           function(x) standardGeneric("anota2seqUtilsGetUTR5"))
setMethod("anota2seqUtilsGetUTR5", "anota2seqUtilsData",
          function(x){
            x@annot@UTR5
          })

setGeneric("anota2seqUtilsGetCDS",
           function(x) standardGeneric("anota2seqUtilsGetCDS"))
setMethod("anota2seqUtilsGetCDS", "anota2seqUtilsData",
          function(x){
            x@annot@CDS
          })

setGeneric("anota2seqUtilsGetID",
           function(x) standardGeneric("anota2seqUtilsGetID"))
setMethod("anota2seqUtilsGetID", "anota2seqUtilsRegion",
          function(x){
            x@id
          })

setGeneric("anota2seqUtilsGetGeneID",
           function(x) standardGeneric("anota2seqUtilsGetGeneID"))
setMethod("anota2seqUtilsGetGeneID", "anota2seqUtilsRegion",
          function(x){
            x@geneID
          })

setGeneric("anota2seqUtilsGetSeq",
           function(x) standardGeneric("anota2seqUtilsGetSeq"))
setMethod("anota2seqUtilsGetSeq", "anota2seqUtilsRegion",
          function(x){
            x@seq
          })

setGeneric("anota2seqUtilsGetDataIn",
           function(x) standardGeneric("anota2seqUtilsGetDataIn"))
setMethod("anota2seqUtilsGetDataIn", "anota2seqUtilsData",
          function(x){
            x@dataIn@geneList
          })

setGeneric("anota2seqUtilsGetBg",
           function(x) standardGeneric("anota2seqUtilsGetBg"))
setMethod("anota2seqUtilsGetBg", "anota2seqUtilsData",
          function(x){
            x@dataIn@background
          })

setGeneric("anota2seqUtilsGetEff",
           function(x) standardGeneric("anota2seqUtilsGetEff"))
setMethod("anota2seqUtilsGetEff", "anota2seqUtilsData",
          function(x){
            x@dataIn@effect
          })

setGeneric("anota2seqUtilsGetColours",
           function(x) standardGeneric("anota2seqUtilsGetColours"))
setMethod("anota2seqUtilsGetColours", "anota2seqUtilsData",
          function(x){
            x@dataIn@colours
          })

setGeneric("anota2seqUtilsGetSpecies",
           function(x) standardGeneric("anota2seqUtilsGetSpecies"))
setMethod("anota2seqUtilsGetSpecies", "anota2seqUtilsData",
          function(x){
            x@species
          })

setGeneric("anota2seqUtilsGetVersion",
           function(x) standardGeneric("anota2seqUtilsGetVersion"))
setMethod("anota2seqUtilsGetVersion", "anota2seqUtilsData",
          function(x){
            x@version
          })

setGeneric("anota2seqUtilsGetSelection",
           function(x) standardGeneric("anota2seqUtilsGetSelection"))
setMethod("anota2seqUtilsGetSelection", "anota2seqUtilsData",
          function(x){
            x@selection
          })

setGeneric("anota2seqUtilsGetGSEA",
           function(x) standardGeneric("anota2seqUtilsGetGSEA"))
setMethod("anota2seqUtilsGetGSEA", "anota2seqUtilsData",
          function(x){
            x@analysis@GSEA
          })

