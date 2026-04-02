mock_region <- function(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("ATG", "CCC")
) {
    methods::new(
        "postNetRegion",
        id = ids,
        geneID = gene_ids,
        sequences = sequences
    )
}

if (!methods::isClass("Anota2seqDataSet")) {
    methods::setClass(
        "Anota2seqDataSet",
        slots = c(
            contrasts = "ANY",
            dataP = "ANY",
            totalmRNA = "ANY",
            translatedmRNA = "ANY",
            buffering = "ANY",
            translation = "ANY"
        )
    )
}

if (!methods::isClass("Anota2seqOutput")) {
    methods::setClass(
        "Anota2seqOutput",
        slots = c(
            apvStatsRvm = "ANY"
        )
    )
}

if (!methods::isClass("mockEnrichResult")) {
    methods::setClass("mockEnrichResult", slots = c(result = "data.frame"))
}

if (!methods::isClass("mockAuc")) {
    methods::setClass("mockAuc", slots = c(y.values = "list"))
}

mock_ads <- function() {
    methods::new(
        "Anota2seqDataSet",
        contrasts = matrix(1, nrow = 1, ncol = 2),
        dataP = matrix(
            1:3,
            ncol = 1,
            dimnames = list(c("gene1", "gene2", "gene3"), "x")
        ),
        totalmRNA = NULL,
        translatedmRNA = NULL,
        buffering = NULL,
        translation = NULL
    )
}

mock_ads_plot <- function() {
    effect_df <- data.frame(apvEff = c(2, -2, 0.5, -0.5), row.names = c("gene1", "gene2", "gene3", "gene4"))
    ads_slots <- methods::getSlots("Anota2seqDataSet")
    data_mat <- matrix(1:4, ncol = 1, dimnames = list(c("gene1", "gene2", "gene3", "gene4"), "x"))
    args <- list(
        Class = "Anota2seqDataSet",
        contrasts = matrix(1, nrow = 1, ncol = 1),
        dataP = data_mat
    )

    if ("totalmRNA" %in% names(ads_slots)) {
        out <- methods::new(
            "Anota2seqOutput",
            apvStatsRvm = list(effect_df)
        )

        args$totalmRNA <- out
        args$translatedmRNA <- out
        args$translation <- out
        args$buffering <- out
    }

    if ("dataT" %in% names(ads_slots)) {
        args$dataT <- data_mat
        args$phenoVec <- c("A", "B", "A", "B")
        args$batchVec <- NULL
        args$qualityControl <- NULL
        args$residOutlierTest <- NULL
        args$selectedTranslatedmRNA <- NULL
        args$selectedTotalmRNA <- NULL
        args$selectedTranslation <- NULL
        args$selectedBuffering <- NULL
        args$mRNAAbundance <- NULL
        args$deltaData <- NULL
        args$regModes <- TRUE
    }

    do.call(methods::new, args)
}

mock_ptn <- function(
        with_background = TRUE,
        with_codon = FALSE,
        with_miRNA = FALSE,
        with_gsea = FALSE,
        with_gage = FALSE,
        with_motifs = FALSE,
        with_feature_models = FALSE
) {
    motif_slot <- if (with_motifs) {
        methods::new(
            "postNetMotifs",
            UTR5 = NULL,
            CDS = list(
                motifSelection = data.frame(score = 1, row.names = "motif1"),
                comparisonA = data.frame(score = 1, row.names = "motif1")
            ),
            UTR3 = list(
                motifSelection = data.frame(score = 2, row.names = "motif2"),
                comparisonA = data.frame(score = 2, row.names = "motif2")
            )
        )
    } else {
        NULL
    }

    codon_slot <- if (with_codon) {
        methods::new(
            "postNetCodons",
            codonAnalysis = methods::new(
                "postNetCodonsAll",
                geneID = c("gene1", "gene1", "gene2"),
                codon = c("AAA", "AAG", "AAA"),
                AA = c("Lys", "Lys", "Lys"),
                count = c(3, 2, 4),
                frequency = c(0.3, 0.2, 0.4),
                AACountPerGene = c(5, 5, 4),
                relative_frequency = c(0.6, 0.4, 1.0)
            ),
            codonSelection = list(comparisonA = c("AAA", "AAG"))
        )
    } else {
        NULL
    }

    mirna_slot <- if (with_miRNA) {
        methods::new(
            "postNetmiRNA",
            miRNA_analysis = list(
                greater = data.frame(
                    score = c(1.1, 1.2),
                    context = c(0.2, 0.4),
                    support = c(11, 12),
                    pvalue = c(0.01, 0.2),
                    qvalue = c(0.02, 0.3),
                    row.names = c("mir1", "mir2")
                ),
                less = data.frame(
                    score = 1.3,
                    context = 0.5,
                    support = 13,
                    pvalue = 0.03,
                    qvalue = 0.04,
                    row.names = "mir3"
                )
            ),
            miRNA_to_gene = list(mir1 = c("gene1", "gene2"), mir3 = "gene2")
        )
    } else {
        NULL
    }

    gsea_slot <- if (with_gsea) {
        data.frame(
            pathway = c("setA", "setB"),
            V2 = 1:2,
            V3 = 3:4,
            V4 = 5:6,
            V5 = 7:8,
            V6 = 9:10,
            V7 = 11:12,
            padj = c(0.01, 0.2)
        )
    } else {
        NULL
    }

    gage_slot <- if (with_gage) {
        methods::new(
            "postNetGAGE",
            BP = list(
                greater = data.frame(
                    stat1 = 1,
                    stat2 = 2,
                    stat3 = 3,
                    stat4 = 4,
                    stat5 = 5,
                    q.val = 0.01,
                    row.names = "termA"
                ),
                less = data.frame(
                    stat1 = 1,
                    stat2 = 2,
                    stat3 = 3,
                    stat4 = 4,
                    stat5 = 5,
                    q.val = 0.2,
                    row.names = "termB"
                )
            ),
            CC = NULL,
            MF = NULL,
            KEGG = NULL
        )
    } else {
        NULL
    }

    feature_models <- if (with_feature_models) {
        list(
            lm = list(methods::new(
                "postNetFeatureIntegration_lm",
                univariateModel = methods::new(
                    "postNetUnivariate",
                    pvalue = c(featureA = 0.01),
                    fdr = c(featureA = 0.02),
                    varianceExplained = c(featureA = 12)
                ),
                stepwiseModel = methods::new(
                    "postNetStepWise",
                    models = list(step1 = list(modelA = "fit")),
                    table = matrix(1, nrow = 1, dimnames = list("featureA", "F"))
                ),
                finalModel = methods::new(
                    "postNetFinalModel",
                    totalVarianceExplained = 25,
                    finalModel = FALSE,
                    table = data.frame(feature = "featureA", value = 1)
                ),
                selectedFeatures = "featureA",
                networkGraph = list(nodes = c("featureA", "effM"))
            )),
            rf = list(methods::new(
                "postNetFeatureIntegration_rf",
                preModel = list(rank = 1),
                borutaModel = list(selected = "featureA"),
                finalModel = list(accuracy = 0.9),
                selectedFeatures = "featureA",
                prediction = c(gene1 = 0.8)
            ))
        )
    } else {
        NULL
    }

    methods::new(
        "postNetData",
        species = "human",
        version = "v1",
        selection = "random",
        seed = 123,
        annot = methods::new(
            "postNetAnnot",
            UTR5 = mock_region(sequences = c("AAA", "GGG")),
            CDS = mock_region(sequences = c("ATGAAA", "ATGAAG")),
            UTR3 = mock_region(sequences = c("TTT", "AAA")),
            CCDS = NULL
        ),
        dataIn = methods::new(
            "postNetDataIn",
            background = if (with_background) c("gene1", "gene2", "gene3") else NULL,
            geneList = list(comparisonA = c("gene1", "gene2")),
            effect = c(gene1 = 1.2, gene2 = -0.4),
            colours = "red"
        ),
        features = data.frame(featureA = c(1, 2), row.names = c("gene1", "gene2")),
        analysis = methods::new(
            "postNetAnalysis",
            featureIntegration = feature_models,
            motifs = motif_slot,
            codons = codon_slot,
            GO = NULL,
            GSEA = gsea_slot,
            GAGE = gage_slot,
            miRNA = mirna_slot
        )
    )
}
