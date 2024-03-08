## anota2seqUtils S4 class implementation
setClassUnion("RegionOrNULL",members=c("anota2seqUtilsRegion", "NULL"))
setClassUnion("characterOrNULL",members=c("character", "NULL"))

setClass("anota2seqUtilsRegion",
         slots = c(
           id = "character",
           geneID = "character",
           seq =  "character"
         )
)

setClass("anota2seqUtilsAnnot",
         slots = c(
           UTR5 = "RegionOrNULL",
           CDS = "RegionOrNULL",
           UTR3 = "RegionOrNULL"
         )
)

setClass("anota2seqUtilsDataIn",
         slots = c(
           background = "characterOrNULL",
           geneList = "list"
         )
)

setClass("anota2seqUtilsFeatures",
         slots = c(
           features = "list"
         )
)

setClass("anota2seqUtilsData",
         slots = c(
           species = "character",
           version = "character",
           annot =  "anota2seqUtilsAnnot",
           dataIn = "anota2seqUtilsDataIn",
           features = "anota2seqUtilsFeatures"
         )
)

