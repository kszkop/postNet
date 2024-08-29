## anota2seqUtils S4 class implementation
setClassUnion("characterOrNULL",members=c("character", "NULL"))
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("dataframeOrNULL",members=c("data.frame", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClassUnion("numericOrNULLOrlogical",members=c("numeric", "NULL", "logical"))
setClassUnion("characterOrNULLOrlogical",members=c("character", "NULL", "logical"))
setClassUnion("anovaOrNULL",members=c("anova", "logical"))

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


setClass("anota2seqUtilsCodonsAll",
         slots = c(
           geneID = "characterOrNULL",
           codon = "characterOrNULL",
           AA = "characterOrNULLOrlogical",
           count = "numericOrNULL",
           frequency = "numericOrNULL",
           AACountPerGene = "numericOrNULLOrlogical",
           relative_frequency = "numericOrNULLOrlogical"
         )
)

setClassUnion("codonsAllOrNULL",members=c("anota2seqUtilsCodonsAll", "NULL"))

setClass("anota2seqUtilsCodons",
         slots = c(
           codonsAll = "codonsAllOrNULL",
           codonsSel = "listOrNULL"
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

setClass("anota2seqUtilsUnivariate",
         slots = c(
           pvalue = "numericOrNULL",
           fdr = "numericOrNULL",
           varianceExplained = "numericOrNULL"
         )
)

setClass("anota2seqUtilsStepWise",
         slots = c(
           models = "listOrNULL",
           table = "matrixOrNULL"
         )
)

setClass("anota2seqUtilsFinalModel",
         slots = c(
           totalVarianceExplained = "numericOrNULL",
           finalModel = "anovaOrNULL",
           table = "dataframeOrNULL"
         )
)

setClassUnion("univariateOrNULL",members=c("anota2seqUtilsUnivariate", "NULL"))
setClassUnion("stepwiseOrNULL",members=c("anota2seqUtilsStepWise", "NULL"))
setClassUnion("finalmodelOrNULL",members=c("anota2seqUtilsFinalModel", "NULL"))

setClass("anota2seqUtilsFeatureIntegration_lm",
         slots = c(
           univariateModel = "univariateOrNULL",
           stepwiseModel = "stepwiseOrNULL",
           finalModel = "finalmodelOrNULL",
           selectedFeatures = "characterOrNULL",
           networkGraph = "ANY"
         )
)

setClass("anota2seqUtilsFeatureIntegration_rf",
         slots = c(
            preModel = "ANY",
            borutaModel = "ANY",
            finalModel = "ANY",
            selectedFeatures = "numericOrNULL"
         )
)


setClassUnion("lmOrNULL",members=c("anota2seqUtilsFeatureIntegration_lm", "NULL"))
setClassUnion("rfOrNULL",members=c("anota2seqUtilsFeatureIntegration_rf", "NULL"))

setClass("anota2seqUtilsFeatureIntegration",
         slots = c(
           lm = "listOrNULL",
           rf  = "listOrNULL"
         )
)

setClassUnion("motifsOrNULL",members=c("anota2seqUtilsMotifs", "NULL"))
setClassUnion("codonsOrNULL",members=c("anota2seqUtilsCodons", "NULL"))
setClassUnion("miRNAOrNULL",members=c("anota2seqUtilsmiRNA", "NULL"))
setClassUnion("GOOrNULL",members=c("anota2seqUtilsGO", "NULL"))
setClassUnion("GAGEOrNULL",members=c("anota2seqUtilsGAGE", "NULL"))
setClassUnion("FIOrNULL",members=c("anota2seqUtilsFeatureIntegration", "NULL"))


setClass("anota2seqUtilsAnalysis",
         slots = c(
           featureIntegration = "FIOrNULL",
           motifs  = "motifsOrNULL",
           codons = "codonsOrNULL",
           GO = "GOOrNULL",
           GSEA = "listOrNULL",
           GAGE = "GAGEOrNULL",
           miRNA = "miRNAOrNULL"
         )
)

setClass("anota2seqUtilsData",
         slots = c(
           species = "character",
           version = "character",
           selection = "character",
           annot =  "anota2seqUtilsAnnot",
           dataIn = "anota2seqUtilsDataIn",
           features = "dataframeOrNULL",
           analysis = "anota2seqUtilsAnalysis"
         )
)

