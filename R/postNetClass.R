## postNet S4 class implementation
setClassUnion("characterOrNULL",members=c("character", "NULL"))
setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("dataframeOrNULL",members=c("data.frame", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))
setClassUnion("numericOrNULLOrlogical",members=c("numeric", "NULL", "logical"))
setClassUnion("characterOrNULLOrlogical",members=c("character", "NULL", "logical"))
setClassUnion("characterOrnumericOrNULL",members=c("character", "numeric", "NULL"))

setClassUnion("anovaOrNULL",members=c("anova", "logical"))

setClass("postNetRegion",
         slots = c(
           id = "character",
           geneID = "character",
           sequences =  "character"
         )
)

setClassUnion("RegionOrNULL",members=c("postNetRegion", "NULL"))

setClass("postNetAnnot",
         slots = c(
           UTR5 = "RegionOrNULL",
           CDS = "RegionOrNULL",
           UTR3 = "RegionOrNULL",
           CCDS = "RegionOrNULL"
         )
)

setClass("postNetDataIn",
         slots = c(
           background = "characterOrNULL",
           geneList = "list",
           effect = "numeric",
           colours = "character"
         )
)

setClass("postNetMotifs",
         slots = c(
           UTR5 = "listOrNULL",
           CDS = "listOrNULL",
           UTR3 = "listOrNULL"
         )
)


setClass("postNetCodonsAll",
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

setClassUnion("codonsAllOrNULL",members=c("postNetCodonsAll", "NULL"))

setClass("postNetCodons",
         slots = c(
           codonAnalysis = "codonsAllOrNULL",
           codonSelection = "listOrNULL"
         )
)

setClass("postNetmiRNA",
         slots = c(
          miRNA_analysis = "listOrNULL",
          miRNA_to_gene = "listOrNULL"
         )
)

setClass("postNetGO",
         slots = c(
           BP = "listOrNULL",
           CC = "listOrNULL",
           MF = "listOrNULL",
           KEGG = "listOrNULL"
         )
)

setClass("postNetGAGE",
         slots = c(
           BP = "listOrNULL",
           CC = "listOrNULL",
           MF = "listOrNULL",
           KEGG = "listOrNULL"
         )
)

setClass("postNetUnivariate",
         slots = c(
           pvalue = "numericOrNULL",
           fdr = "numericOrNULL",
           varianceExplained = "numericOrNULL"
         )
)

setClass("postNetStepWise",
         slots = c(
           models = "listOrNULL",
           table = "matrixOrNULL"
         )
)

setClass("postNetFinalModel",
         slots = c(
           totalVarianceExplained = "numericOrNULL",
           finalModel = "anovaOrNULL",
           table = "dataframeOrNULL"
         )
)

setClassUnion("univariateOrNULL",members=c("postNetUnivariate", "NULL"))
setClassUnion("stepwiseOrNULL",members=c("postNetStepWise", "NULL"))
setClassUnion("finalmodelOrNULL",members=c("postNetFinalModel", "NULL"))

setClass("postNetFeatureIntegration_lm",
         slots = c(
           univariateModel = "univariateOrNULL",
           stepwiseModel = "stepwiseOrNULL",
           finalModel = "finalmodelOrNULL",
           selectedFeatures = "characterOrnumericOrNULL",
           networkGraph = "ANY"
         )
)

setClass("postNetFeatureIntegration_rf",
         slots = c(
            preModel = "ANY",
            borutaModel = "ANY",
            finalModel = "ANY",
            selectedFeatures = "characterOrnumericOrNULL",
            prediction = "ANY"
         )
)


#setClassUnion("lmOrNULL",members=c("postNetFeatureIntegration_lm", "NULL"))
#setClassUnion("rfOrNULL",members=c("postNetFeatureIntegration_rf", "NULL"))

#setClass("postNetFeatureIntegration",
#         slots = c(
#           lm = "listOrNULL",
#           rf  = "listOrNULL",
#           featureMap = "ANY"
#         )
#)

setClassUnion("motifsOrNULL",members=c("postNetMotifs", "NULL"))
setClassUnion("codonsOrNULL",members=c("postNetCodons", "NULL"))
setClassUnion("miRNAOrNULL",members=c("postNetmiRNA", "NULL"))
setClassUnion("GOOrNULL",members=c("postNetGO", "NULL"))
setClassUnion("GAGEOrNULL",members=c("postNetGAGE", "NULL"))
setClassUnion("FIOrNULL",members=c("postNetFeatureIntegration", "NULL"))


setClass("postNetAnalysis",
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

setClass("postNetData",
         slots = c(
           species = "characterOrNULL",
           version = "characterOrNULL",
           selection = "character",
           seed = "ANY",
           annot =  "postNetAnnot",
           dataIn = "postNetDataIn",
           features = "dataframeOrNULL",
           analysis = "postNetAnalysis"
         )
)

