## anota2seqUtils S4 class implementation
setClassUnion("RegionOrNULL",members=c("anota2seqUtilsRegion", "NULL"))
setClassUnion("miRNAOrNULL",members=c("anota2seqUtilsmiRNA", "NULL"))
setClassUnion("GOOrNULL",members=c("anota2seqUtilsGO", "NULL"))
setClassUnion("GAGEOrNULL",members=c("anota2seqUtilsGAGE", "NULL"))
setClassUnion("characterOrNULL",members=c("character", "NULL"))
setClassUnion("listOrNULL",members=c("list", "NULL"))

setClass("anota2seqUtilsRegion",
         slots = c(
           id = "characterOrNULL",
           geneID = "characterOrNULL",
           seq =  "characterOrNULL"
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
           geneList = "list",
           effect = "numeric",
           colours = "character"
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

setClass("anota2seqUtilsAnalysis",
         slots = c(
           featureIntegration = "listOrNULL",
           motifs  = "listOrNULL",
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

