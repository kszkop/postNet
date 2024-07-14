## anota2seqUtils S4 class implementation
setClassUnion("characterOrNULL",members=c("character", "NULL"))
setClassUnion("listOrNULL",members=c("list", "NULL"))

setClass("anota2seqUtilsRegion",
         slots = c(
           id = "character",
           geneID = "character",
           seq =  "character"
         )
)

setClassUnion("RegionOrNULL",members=c("anota2seqUtilsRegion", "NULL"))

setClass("anota2seqUtilsAnnot",
         slots = c(
           UTR5 = "RegionOrNULL",
           CDS = "RegionOrNULL",
           UTR3 = "RegionOrNULL",
           CCDS = "RegionOrNULL"
         )
)

setClass("anota2seqUtilsDataIn",
         slots = c(
           background = "characterOrNULL",
           geneList = "list",
           effect = "numeric",
           colours = "character"
         )
)

setClass("anota2seqUtilsMotifs",
         slots = c(
           UTR5 = "listOrNULL",
           CDS = "listOrNULL",
           UTR3 = "listOrNULL"
         )
)

setClass("anota2seqUtilsmiRNA",
         slots = c(
          miRNA_analysis = "listOrNULL",
          miRNA_to_gene = "listOrNULL"
         )
)

setClass("anota2seqUtilsGO",
         slots = c(
           BP = "listOrNULL",
           CC = "listOrNULL",
           MF = "listOrNULL",
           KEGG = "listOrNULL"
         )
)

setClass("anota2seqUtilsGAGE",
         slots = c(
           BP = "listOrNULL",
           CC = "listOrNULL",
           MF = "listOrNULL",
           KEGG = "listOrNULL"
         )
)

setClassUnion("motifsOrNULL",members=c("anota2seqUtilsMotifs", "NULL"))
setClassUnion("miRNAOrNULL",members=c("anota2seqUtilsmiRNA", "NULL"))
setClassUnion("GOOrNULL",members=c("anota2seqUtilsGO", "NULL"))
setClassUnion("GAGEOrNULL",members=c("anota2seqUtilsGAGE", "NULL"))

setClass("anota2seqUtilsAnalysis",
         slots = c(
           featureIntegration = "listOrNULL",
           motifs  = "motifsOrNULL",
           codons = "listOrNULL",
           GO = "GOOrNULL",
           GSEA = "listOrNULL",
           GAGE = "GAGEOrNULL",
           miRNA = "miRNAOrNULL"
         )
)

setClass("anota2seqUtilsFeatures",
         slots = c(
           features = "listOrNULL"
         )
)

setClass("anota2seqUtilsData",
         slots = c(
           species = "character",
           version = "character",
           selection = "character",
           annot =  "anota2seqUtilsAnnot",
           dataIn = "anota2seqUtilsDataIn",
           features = "anota2seqUtilsFeatures",
           analysis = "anota2seqUtilsAnalysis"
         )
)

