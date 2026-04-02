test_that("lengthAnalysis and contentAnalysis return deterministic summaries", {
    ptn <- mock_ptn(with_background = TRUE)

    length_res <- lengthAnalysis(ptn, region = "UTR3", plotOut = FALSE)
    expect_named(length_res, "UTR3_length")
    expect_equal(unname(length_res$UTR3_length), rep(log2(3), 2))
    expect_equal(names(length_res$UTR3_length), c("gene1", "gene2"))

    content_res <- contentAnalysis(
        ptn,
        contentIn = c("A"),
        region = "UTR3",
        plotOut = FALSE
    )
    expect_named(content_res, "UTR3_A")
    expect_equal(names(content_res$UTR3_A), c("gene1", "gene2"))
    expect_equal(round(unname(content_res$UTR3_A), 3), c(0, 100))

    pos_res <- contentAnalysis(
        ptn,
        contentIn = c("A1"),
        region = "CDS",
        plotOut = FALSE
    )
    expect_named(pos_res, "CDS_A1")
    expect_equal(round(unname(pos_res$CDS_A1), 3), c(33.333, 33.333))
})

test_that("contentAnalysis and lengthAnalysis validate plotting and comparison inputs", {
    ptn_no_bg <- mock_ptn(with_background = FALSE)

    expect_error(
        lengthAnalysis(ptn_no_bg, region = "UTR3", comparisons = list(c(0, 1)), plotOut = FALSE),
        "0 always denotes the background"
    )
    expect_error(
        contentAnalysis(ptn_no_bg, contentIn = "A", region = "UTR3", comparisons = list("bad"), plotOut = FALSE),
        "must be a list of numeric vectors"
    )
    expect_error(
        contentAnalysis(ptn_no_bg, contentIn = "A", region = "UTR3", subregion = 2, plotOut = FALSE),
        "Please specify if you would like to 'select' or 'exclude'"
    )
})

test_that("lengthAnalysis, contentAnalysis and codonCalc write plot outputs", {
    ptn <- mock_ptn(with_background = TRUE, with_codon = TRUE)
    out_dir <- tempfile("plots-basic")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    testthat::local_mocked_bindings(
        plotPostNet = function(...) NULL,
        .package = "postNet"
    )

    len_plot <- lengthAnalysis(
        ptn,
        region = "UTR3",
        comparisons = list(c(0, 1)),
        plotOut = TRUE,
        pdfName = "lenplot"
    )
    expect_named(len_plot, "UTR3_length")
    expect_true(file.exists(file.path(out_dir, "lenplot_UTR3_boxplot_lengthAnalysis.pdf")))

    content_plot <- contentAnalysis(
        ptn,
        contentIn = "A",
        region = "UTR3",
        comparisons = list(c(0, 1)),
        plotOut = TRUE,
        pdfName = "contentplot"
    )
    expect_named(content_plot, "UTR3_A")
    expect_true(file.exists(file.path(out_dir, "contentplot_UTR3_A_content.pdf")))

    codon_plot <- codonCalc(
        ptn,
        featsel = list(lys = "K"),
        analysis = "AA",
        unit = "freq",
        comparisons = list(c(0, 1)),
        plotOut = TRUE,
        pdfName = "codoncalc"
    )
    expect_named(codon_plot, "lys")
    expect_true(file.exists(file.path(out_dir, "codoncalc_features_lys_codonCalc.pdf")))
})

test_that("contentAnalysis covers subregion handling and constant-value plotting path", {
    ptn <- mock_ptn(with_background = TRUE)
    ptn@dataIn@background <- c("gene1", "gene2")
    methods::slot(ptn@annot, "UTR3") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("AAAA", "AAAA")
    )

    sub_sel <- contentAnalysis(
        ptn,
        contentIn = "A",
        region = "UTR3",
        subregion = 2,
        subregionSel = "select",
        plotOut = FALSE
    )
    expect_equal(sub_sel$UTR3_A, c(gene1 = 100, gene2 = 100))

    out_dir <- tempfile("content-const")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)
    expect_message(
        contentAnalysis(
            ptn,
            contentIn = "A",
            region = "UTR3",
            comparisons = list(c(0, 1)),
            plotOut = TRUE,
            pdfName = "const"
        ),
        "No plot will be produced as all values are the same"
    )
    expect_false(any(grepl("^const_.*content\\.pdf$", list.files(out_dir))))
})

test_that("codonCalc summarizes codon counts and frequencies without plotting", {
    ptn <- mock_ptn(with_background = TRUE, with_codon = TRUE)

    count_res <- codonCalc(
        ptn,
        featsel = list(lys = c("AAA", "AAG")),
        analysis = "codon",
        unit = "count",
        plotOut = FALSE
    )
    expect_equal(count_res$lys, c(gene1 = 5, gene2 = 4))

    freq_res <- codonCalc(
        ptn,
        featsel = list(lys = c("AAA", "AAG")),
        analysis = "codon",
        unit = "freq",
        plotOut = FALSE
    )
    expect_equal(freq_res$lys, c(gene1 = 0.5, gene2 = 0.4))

    expect_error(
        codonCalc(ptn, featsel = list("AAA"), plotOut = FALSE),
        "must be a named list"
    )
    expect_error(
        codonCalc(ptn, featsel = list(lys = "AAA"), unit = "counts", plotOut = FALSE),
        "must be either 'count' or 'freq'"
    )
    expect_error(
        codonCalc(mock_ptn(), featsel = list(lys = "AAA"), plotOut = FALSE),
        "no applicable method"
    )
    expect_error(
        codonCalc(ptn, featsel = list(lys = "AAA"), plotOut = "yes"),
        "plotOut"
    )
    expect_error(
        codonCalc(ptn, featsel = list(lys = "AAA"), plotOut = TRUE, plotType = NULL),
        "Please provide a selection for 'plotType'"
    )
    expect_error(
        codonCalc(
            mock_ptn(with_background = FALSE, with_codon = TRUE),
            featsel = list(lys = "AAA"),
            comparisons = list(c(0, 1)),
            plotOut = FALSE
        ),
        "0 always denotes the background"
    )
    expect_error(
        codonCalc(ptn, featsel = list(lys = "AAA"), analysis = "protein", plotOut = FALSE),
        "The options are: 'codon' or 'AA'"
    )
    expect_error(
        codonCalc(ptn, featsel = list(lys = "XYZ"), analysis = "codon", plotOut = FALSE),
        "Invalid codons provided"
    )
    aa_res <- codonCalc(
        ptn,
        featsel = list(lys = c("K", "Lys")),
        analysis = "AA",
        unit = "count",
        plotOut = FALSE
    )
    expect_equal(aa_res$lys, c(gene1 = 5, gene2 = 4))
})

test_that("feature integration accessors expose stored models and selections", {
    ptn <- mock_ptn(with_feature_models = TRUE)

    expect_equal(ptn_selectedFeatures(ptn, analysis_type = "lm", comparison = 1), "featureA")
    expect_equal(ptn_selectedFeatures(ptn, analysis_type = "rf", comparison = 1), "featureA")
    expect_equal(ptn_networkGraph(ptn, comparison = 1), list(nodes = c("featureA", "effM")))

    final_model <- ptn_model(
        ptn,
        analysis_type = "lm",
        model = "finalModel",
        comparison = 1
    )
    expect_s4_class(final_model, "postNetFinalModel")

    rf_model <- ptn_model(
        ptn,
        analysis_type = "rf",
        model = "finalModel",
        comparison = 1
    )
    expect_equal(rf_model, list(accuracy = 0.9))

    expect_error(
        ptn_model(ptn, analysis_type = "lm", model = "finalModel", comparison = 2),
        "There are only 1 comparisons"
    )
})

test_that("slopeFilt validates inputs and filters mocked slopes", {
    ads <- mock_ads()

    expect_error(
        slopeFilt(ads, regulationGen = "abundance", contrastSel = 1),
        "should be either 'translation' or 'buffering'"
    )
    expect_error(
        slopeFilt(ads, regulationGen = "translation", contrastSel = 3),
        "should be a number corresponding to the desired contrast"
    )

    testthat::local_mocked_bindings(
        anota2seqGetOutput = function(...) {
            data.frame(
                slope = c(-2, 0.5, 3),
                other = c(1, 1, 1),
                row.names = c("gene1", "gene2", "gene3")
            )
        },
        .package = "anota2seq"
    )

    expect_equal(
        slopeFilt(ads, regulationGen = "translation", contrastSel = 1),
        c("gene1", "gene3")
    )
    expect_equal(
        slopeFilt(
            ads,
            regulationGen = "buffering",
            contrastSel = 1,
            minSlope = -3,
            maxSlope = 2
        ),
        "gene3"
    )
})

test_that("miRNAanalysis validates thresholds and stores mocked enrichment", {
    ptn <- mock_ptn(with_background = TRUE)
    targetscan <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            Cumulative.weighted.context...score = c(-0.4, -0.3),
            Aggregate.PCT = c(0.4, 0.5),
            Gene.Symbol = c("gene1", "gene2"),
            Representative.miRNA = c("mirA", "mirB"),
            Species.ID = c(9606, 9606)
        ),
        file = targetscan,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    expect_error(
        miRNAanalysis(ptn, targetscan, maxSize = 5, minSize = 5),
        "maxSize must be greater than minSize"
    )
    expect_error(
        miRNAanalysis(ptn, targetscan, contextScore = 0.1),
        "negative value for 'contextScore'"
    )
    expect_error(
        miRNAanalysis(ptn, targetscan, Pct = 2),
        "value between 0 and 1 for 'Pct'"
    )

    testthat::local_mocked_bindings(
        gage = function(exprs, gsets, set.size, rank.test, use.fold) {
            expect_equal(sort(names(gsets)), c("mirA", "mirB"))
            list(
                greater = data.frame(
                    stat.mean = 1,
                    p.val = 0.01,
                    set.size = 1,
                    q.val = 0.02,
                    rank = 1,
                    row.names = "mirA"
                ),
                less = data.frame(
                    stat.mean = 0.5,
                    p.val = 0.2,
                    set.size = 1,
                    q.val = 0.3,
                    rank = 1,
                    row.names = "mirB"
                )
            )
        },
        .package = "gage"
    )

    out <- miRNAanalysis(
        ptn,
        targetscan,
        genesSlopeFiltOut = "gene2",
        contextScore = -0.2,
        Pct = 0.1,
        maxSize = 50,
        minSize = 1
    )

    expect_s4_class(out@analysis@miRNA, "postNetmiRNA")
    expect_equal(out@analysis@miRNA@miRNA_to_gene$mirA, "gene1")
    expect_equal(rownames(out@analysis@miRNA@miRNA_analysis$greater), "mirA")
})

test_that("motifAnalysis validates setup and stores mocked STREME output", {
    ptn_no_bg <- mock_ptn(with_background = FALSE)
    expect_error(
        motifAnalysis(ptn_no_bg, memePath = "/tmp/meme", region = "UTR3"),
        "background gene set must be provided"
    )

    ptn <- mock_ptn(with_background = TRUE)
    expect_error(
        motifAnalysis(ptn, memePath = NULL, region = "UTR3"),
        "Please provide full file path to the STREME executables"
    )
    expect_error(
        motifAnalysis(ptn, memePath = "/tmp/meme", seqType = "lipid", region = "UTR3"),
        "must be either 'dna', 'rna', or 'protein'"
    )

    testthat::local_mocked_bindings(
        runStreme = function(input, control, meme_path, alph, minw) {
            data.frame(
                consensus = c("AUGCUA", "UUUGGG"),
                pval = c(0.001, 0.2)
            )
        },
        .package = "memes"
    )

    out <- motifAnalysis(
        ptn,
        stremeThreshold = 0.05,
        minwidth = 6,
        memePath = "/tmp/meme",
        seqType = "dna",
        region = "UTR3"
    )

    expect_s4_class(out@analysis@motifs, "postNetMotifs")
    expect_equal(out@analysis@motifs@UTR3$motifSelection, "AUGCUA")
    expect_equal(rownames(out@analysis@motifs@UTR3$comparisonA), "1")
})

test_that("motifAnalysis covers empty STREME and subregion filtering branches", {
    ptn <- mock_ptn(with_background = TRUE)

    testthat::local_mocked_bindings(
        runStreme = function(input, control, meme_path, alph, minw) {
            data.frame(consensus = character(), pval = numeric())
        },
        .package = "memes"
    )

    expect_message(
        out <- motifAnalysis(
            ptn,
            stremeThreshold = 0.05,
            minwidth = 6,
            memePath = "/tmp/meme",
            seqType = "dna",
            region = "UTR3",
            subregion = 10,
            subregionSel = "select"
        ),
        "selected subregion is longer than the reference sequence region"
    )
    expect_s4_class(out@analysis@motifs, "postNetMotifs")
    expect_equal(length(out@analysis@motifs@UTR3$motifSelection), 0)
})

test_that("gseaAnalysis and gageAnalysis work with mocked enrichment backends", {
    ptn <- mock_ptn(with_background = TRUE)

    expect_error(
        gseaAnalysis(ptn, geneSet = NULL, collection = NULL),
        "either 'geneSet' or 'collection'"
    )
    expect_error(
        gseaAnalysis(ptn, geneSet = list(setA = "gene1"), maxSize = 1, minSize = 1),
        "must be greater than 'minSize'"
    )

    testthat::local_mocked_bindings(
        fgsea = function(pathways, stats, minSize, maxSize) {
            data.frame(
                pathway = "setA",
                pval = 0.01,
                padj = 0.02,
                log2err = 0.1,
                ES = 1.5,
                NES = 1.2,
                size = 2,
                leadingEdge = I(list(c("gene1", "gene2")))
            )
        },
        .package = "fgsea"
    )

    out_gsea <- gseaAnalysis(
        ptn,
        geneSet = list(setA = c("gene1", "gene2")),
        maxSize = 10,
        minSize = 1,
        name = tempfile("gsea")
    )

    expect_equal(out_gsea@analysis@GSEA$Term, "setA")
    expect_equal(out_gsea@analysis@GSEA$Genes, "gene1:gene2")
    expect_error(
        gseaAnalysis(
            ptn,
            genesSlopeFiltOut = names(ptn_effect(ptn)),
            geneSet = list(setA = c("gene1", "gene2")),
            maxSize = 10,
            minSize = 1
        ),
        "No ranked genes are available"
    )

    testthat::local_mocked_bindings(
        ExperimentHub = function(...) structure(list(), class = "mock_eh"),
        .package = "ExperimentHub"
    )
    testthat::local_mocked_bindings(
        query = function(x, pattern) x,
        .package = "AnnotationHub"
    )
    testthat::local_mocked_bindings(
        getMsigdbVersions = function() c("7.5", "7.4"),
        getMsigdb = function(org, id, version) list(setA = c("gene1", "gene2")),
        appendKEGG = function(x, version) x,
        subsetCollection = function(x, collection, subcollection = NULL) x,
        .package = "msigdb"
    )
    testthat::local_mocked_bindings(
        geneIds = function(x) x,
        .package = "GSEABase"
    )
    testthat::local_mocked_bindings(
        fgsea = function(pathways, stats, minSize, maxSize) {
            data.frame(
                pathway = character(),
                pval = numeric(),
                padj = numeric(),
                log2err = numeric(),
                ES = numeric(),
                NES = numeric(),
                size = integer(),
                leadingEdge = I(list())
            )
        },
        .package = "fgsea"
    )
    out_empty <- gseaAnalysis(
        ptn,
        collection = "h",
        maxSize = 10,
        minSize = 1,
        name = tempfile("gseaempty")
    )
    expect_true(is.data.frame(out_empty@analysis@GSEA))
    expect_equal(nrow(out_empty@analysis@GSEA), 0)

    testthat::local_mocked_bindings(
        convertSymbolToEntrezID = function(geneList, species) {
            stats::setNames(as.character(seq_along(geneList)), geneList)
        },
        convertEntrezIDToSymbol = function(entrezIDList, species) {
            paste0("gene", entrezIDList)
        },
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        go.gsets = function(species, id.type) {
            list(
                go.sets = list(pathA = c("1", "2"), pathB = "2"),
                go.subs = list(BP = c("pathA", "pathB"), MF = character(), CC = character())
            )
        },
        gage = function(exprs, gsets, set.size, rank.test, use.fold) {
            list(
                greater = data.frame(
                    stat.mean = 1,
                    p.val = 0.01,
                    set.size = 2,
                    q.val = 0.02,
                    rank = 1,
                    row.names = "pathA"
                ),
                less = data.frame()
            )
        },
        .package = "gage"
    )

    out_gage <- gageAnalysis(
        ptn,
        category = "BP",
        maxSize = 10,
        minSize = 1
    )

    expect_s4_class(out_gage@analysis@GAGE, "postNetGAGE")
    expect_equal(out_gage@analysis@GAGE@BP$greater$id, "pathA")
    expect_equal(out_gage@analysis@GAGE@BP$greater$Genes, "gene1:gene2")
})

test_that("gageAnalysis covers validation, KEGG and less-result branches", {
    ptn <- mock_ptn(with_background = TRUE)

    expect_error(gageAnalysis(ptn, category = "BP", maxSize = 0, minSize = 1), "must be positive")
    expect_error(gageAnalysis(ptn, category = "BP", maxSize = "big", minSize = 1), "must be integers")

    ptn_mouse <- ptn
    ptn_mouse@species <- "mouse"
    testthat::local_mocked_bindings(
        convertSymbolToEntrezID = function(geneList, species) {
            stats::setNames(as.character(seq_along(geneList)), geneList)
        },
        convertEntrezIDToSymbol = function(entrezIDList, species) {
            paste0("gene", entrezIDList)
        },
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        go.gsets = function(species, id.type) {
            list(
                go.sets = list(pathBP = c("1", "2"), pathMF = c("2", "3"), pathCC = c("3")),
                go.subs = list(BP = "pathBP", MF = "pathMF", CC = "pathCC")
            )
        },
        kegg.gsets = function(species, id.type) {
            list(kg.sets = list(pathKEGG = c("1", "3")))
        },
        gage = function(exprs, gsets, set.size, rank.test, use.fold) {
            list(
                greater = data.frame(
                    stat.mean = 1,
                    p.val = 0.01,
                    set.size = 2,
                    q.val = 0.02,
                    rank = 1,
                    row.names = names(gsets)[1]
                ),
                less = data.frame(
                    stat.mean = -1,
                    p.val = 0.03,
                    set.size = 1,
                    q.val = 0.04,
                    rank = 2,
                    row.names = names(gsets)[1]
                )
            )
        },
        .package = "gage"
    )

    out_multi <- gageAnalysis(
        ptn_mouse,
        category = c("BP", "MF", "CC", "KEGG"),
        genesSlopeFiltOut = "gene2",
        maxSize = 10,
        minSize = 1
    )
    expect_equal(out_multi@analysis@GAGE@BP$greater$Genes, "gene1")
    expect_equal(out_multi@analysis@GAGE@MF$less$id, "pathMF")
    expect_equal(out_multi@analysis@GAGE@CC$greater$id, "pathCC")
    expect_equal(out_multi@analysis@GAGE@KEGG$less$Genes, "gene1")
})

test_that("foldingEnergyAnalysis supports custom folding energy inputs", {
    ptn <- mock_ptn(with_background = TRUE)
    custom_fe <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            id = c("tx1", "tx2"),
            fold_energy = c(-10, -20),
            length = c(6, 6)
        ),
        file = custom_fe,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    expect_error(
        foldingEnergyAnalysis(ptn, sourceFE = "custom", customFileFE = NULL, region = "UTR3", plotOut = FALSE),
        "Please provide a custom file"
    )
    expect_error(
        foldingEnergyAnalysis(ptn, sourceFE = "custom", customFileFE = custom_fe, residFE = "yes", region = "UTR3", plotOut = FALSE),
        "must be logical"
    )
    expect_error(
        foldingEnergyAnalysis(ptn, sourceFE = "remote", region = "UTR3", plotOut = FALSE),
        "Invalid 'sourceFE'"
    )
    expect_error(
        foldingEnergyAnalysis(ptn, sourceFE = "custom", customFileFE = custom_fe, region = "UTR3", plotOut = "yes"),
        "plotOut"
    )
    expect_error(
        foldingEnergyAnalysis(ptn, sourceFE = "custom", customFileFE = custom_fe, region = "UTR3", plotOut = TRUE, plotType = NULL),
        "Please provide an input for 'plotType'"
    )

    fe_res <- foldingEnergyAnalysis(
        ptn,
        sourceFE = "custom",
        customFileFE = custom_fe,
        residFE = FALSE,
        plotOut = FALSE,
        region = "UTR3"
    )
    expect_equal(fe_res$custom, c(gene1 = -10, gene2 = -20))

    fe_resid <- foldingEnergyAnalysis(
        ptn,
        sourceFE = "custom",
        customFileFE = custom_fe,
        residFE = TRUE,
        plotOut = FALSE,
        region = "UTR3"
    )
    expect_equal(round(sum(fe_resid$custom), 8), 0)

    testthat::local_mocked_bindings(
        read.delim = function(...) data.frame(id = c("tx1", "tx2"), fold_energy = c(-5, -7), length = c(6, 6)),
        .package = "utils"
    )
    ptn@species <- "human"
    expect_error(
        foldingEnergyAnalysis(ptn, sourceFE = "load", residFE = FALSE, plotOut = FALSE),
        "Please provide an input for 'region'"
    )
    fe_load <- foldingEnergyAnalysis(
        ptn,
        sourceFE = "load",
        region = "UTR3",
        residFE = FALSE,
        plotOut = FALSE
    )
    expect_named(fe_load, "UTR3_foldingEnergy")

    out_dir <- tempfile("feplots")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)
    testthat::local_mocked_bindings(
        plotPostNet = function(...) NULL,
        .package = "postNet"
    )
    fe_plot <- foldingEnergyAnalysis(
        ptn,
        sourceFE = "custom",
        customFileFE = custom_fe,
        residFE = FALSE,
        region = "UTR3",
        comparisons = list(c(0, 1)),
        plotOut = TRUE,
        plotType = "boxplot",
        pdfName = "feplot"
    )
    expect_named(fe_plot, "custom")
    expect_true(file.exists(file.path(out_dir, "feplot_custom_boxplot_foldEnergyAnalysis.pdf")))

    fe_load_plot <- foldingEnergyAnalysis(
        ptn,
        sourceFE = "load",
        region = "UTR3",
        residFE = FALSE,
        comparisons = list(c(0, 1)),
        plotOut = TRUE,
        plotType = "ecdf",
        pdfName = "feload"
    )
    expect_named(fe_load_plot, "UTR3_foldingEnergy")
    expect_true(file.exists(file.path(out_dir, "feload_UTR3_ecdf_foldEnergyAnalysis.pdf")))

    ptn_bad_species <- ptn
    ptn_bad_species@species <- "rat"
    expect_error(
        foldingEnergyAnalysis(ptn_bad_species, sourceFE = "load", region = "UTR3", plotOut = FALSE),
        "not supported"
    )
})

test_that("contentMotifs validates inputs and counts motif occurrences", {
    ptn <- mock_ptn(with_background = TRUE)
    ptn_no_bg <- mock_ptn(with_background = FALSE)

    expect_error(
        contentMotifs(ptn_no_bg, motifsIn = "AAA", region = "UTR3", plotOut = FALSE, comparisons = list(c(0, 1))),
        "0 always denotes the background"
    )
    expect_error(
        contentMotifs(ptn, motifsIn = "AAA", region = "UTR3", plotOut = FALSE, seqType = "lipid"),
        "must be either 'dna', 'rna', or 'protein'"
    )
    expect_error(
        contentMotifs(ptn, motifsIn = 1, region = "UTR3", plotOut = FALSE),
        "should be a character vector of sequence motifs"
    )

    motif_counts <- contentMotifs(
        ptn,
        motifsIn = "A",
        region = "UTR3",
        unitOut = "number",
        plotOut = FALSE
    )
    expect_named(motif_counts, "UTR3_A")
    expect_equal(motif_counts$UTR3_A, c(gene1 = 0, gene2 = 1))

    motif_pos <- contentMotifs(
        ptn,
        motifsIn = "A",
        region = "UTR3",
        unitOut = "position",
        plotOut = FALSE
    )
    expect_true(all(is.na(motif_pos$UTR3_A$gene1)))
    expect_equal(motif_pos$UTR3_A$gene2$start, 1)
    expect_equal(motif_pos$UTR3_A$gene2$end, 3)
})

test_that("contentMotifs covers protein, G4 and plotting branches", {
    ptn <- mock_ptn(with_background = TRUE)
    methods::slot(ptn@annot, "CDS") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("ATGAAA", "ATGAAG")
    )

    protein_counts <- contentMotifs(
        ptn,
        motifsIn = "M",
        seqType = "protein",
        region = "CDS",
        unitOut = "number",
        plotOut = FALSE
    )
    expect_equal(protein_counts$CDS_M, c(gene1 = 1, gene2 = 1))

    testthat::local_mocked_bindings(
        calc_g4 = function(seq, min_score, unit) 2,
        plotPostNet = function(...) NULL,
        .package = "postNet"
    )
    out_dir <- tempfile("motifplots")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    g4_counts <- contentMotifs(
        ptn,
        motifsIn = "G4",
        region = "UTR3",
        unitOut = "number",
        plotOut = TRUE,
        comparisons = list(c(0, 1)),
        pdfName = "g4plot"
    )
    expect_equal(g4_counts$UTR3_G4, c(gene1 = 2, gene2 = 2))
    expect_true(file.exists(file.path(out_dir, "g4plot_UTR3_G4_content.pdf")))
})

test_that("contentMotifs covers additional validation and edge branches", {
    ptn <- mock_ptn(with_background = TRUE)

    expect_error(
        contentMotifs(ptn, motifsIn = "A", region = "UTR3", plotOut = "yes"),
        "must be logical"
    )
    expect_error(
        contentMotifs(ptn, motifsIn = "A", region = "UTR3", plotOut = FALSE, comparisons = "bad"),
        "must be a list of numeric"
    )
    expect_error(
        contentMotifs(ptn, motifsIn = "A", region = "UTR3", plotOut = FALSE, unitOut = "freq"),
        "must be either 'number' or 'position'"
    )
    expect_error(
        contentMotifs(ptn, motifsIn = "A", region = "UTR3", plotOut = FALSE, resid = "yes"),
        "must be logical"
    )
    expect_error(
        contentMotifs(ptn, motifsIn = "G4", region = "UTR3", plotOut = FALSE, min_score = "high"),
        "minimal score"
    )

    ptn_bad_cds <- ptn
    methods::slot(ptn_bad_cds@annot, "CDS") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("ATGAA", "ATGAAG")
    )
    expect_error(
        contentMotifs(
            ptn_bad_cds,
            motifsIn = "M",
            seqType = "protein",
            region = "CDS",
            unitOut = "number",
            plotOut = FALSE
        ),
        "divided into codons"
    )

    expect_message(
        no_sites <- contentMotifs(
            ptn,
            motifsIn = "CCC",
            region = "UTR3",
            unitOut = "position",
            plotOut = FALSE
        ),
        "does not have any sites"
    )
    expect_equal(no_sites, list())

    expect_message(
        removed <- contentMotifs(
            ptn,
            motifsIn = "A",
            region = "UTR3",
            unitOut = "number",
            subregion = 10,
            subregionSel = "select",
            plotOut = FALSE
        ),
        "selected subregion is longer"
    )
    expect_equal(removed, list())

    excluded <- contentMotifs(
        ptn,
        motifsIn = "A",
        region = "UTR5",
        unitOut = "number",
        subregion = 2,
        subregionSel = "exclude",
        plotOut = FALSE
    )
    expect_named(excluded, "UTR5_A")
    expect_true(excluded$UTR5_A[["gene1"]] >= excluded$UTR5_A[["gene2"]])
})

test_that("goAnalysis validates setup and reshapes mocked enrichment results", {
    ptn <- mock_ptn(with_background = TRUE)

    expect_error(
        goAnalysis(ptn, category = "BP", maxSize = 1, minSize = 1),
        "must be greater than 'minSize'"
    )

    ptn_bad_species <- ptn
    ptn_bad_species@species <- "rat"
    expect_error(
        goAnalysis(ptn_bad_species, category = "BP"),
        "currently only available for human and mouse"
    )

    testthat::local_mocked_bindings(
        convertSymbolToEntrezID = function(geneList, species) {
            stats::setNames(as.character(seq_along(geneList)), geneList)
        },
        convertEntrezIDToSymbol = function(entrezIDList, species) {
            paste0("gene", entrezIDList)
        },
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        enrichGO = function(...) {
            methods::new(
                "mockEnrichResult",
                result = data.frame(
                    ID = "GO:0001",
                    Description = "translation",
                    Count = 5,
                    BgRatio = "20/100",
                    pvalue = 0.001,
                    geneID = "1/2/3/4/5"
                )
            )
        },
        enrichKEGG = function(...) {
            methods::new(
                "mockEnrichResult",
                result = data.frame(
                    ID = "mmu0001",
                    Description = "kegg",
                    category = "Metabolism",
                    subcategory = "Global",
                    Count = 5,
                    BgRatio = "20/100",
                    pvalue = 0.001,
                    geneID = "1/2/3/4/5"
                )
            )
        },
        .package = "clusterProfiler"
    )
    testthat::local_mocked_bindings(
        WriteXLS = function(...) NULL,
        .package = "WriteXLS"
    )

    out <- goAnalysis(
        ptn,
        category = c("BP", "KEGG"),
        maxSize = 10,
        minSize = 1,
        counts = 1,
        FDR = 0.05,
        name = tempfile("go")
    )

    expect_s4_class(out@analysis@GO, "postNetGO")
    expect_equal(out@analysis@GO@BP$comparisonA@result$geneID, "gene1:gene2:gene3:gene4:gene5")
    expect_equal(out@analysis@GO@BP$comparisonA@result$Size, 20)
    expect_equal(out@analysis@GO@KEGG$comparisonA@result$category, "Metabolism")
})

test_that("uorfAnalysis validates inputs and returns deterministic counts", {
    ptn <- mock_ptn(with_background = TRUE)
    methods::slot(ptn@annot, "UTR5") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("AAAATGG", "CCCCCCC")
    )
    methods::slot(ptn@annot, "CDS") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("TAA", "GGG")
    )
    methods::slot(ptn@annot, "UTR3") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("AAA", "TTT")
    )

    expect_error(
        uorfAnalysis(ptn, startCodon = "AT", plotOut = FALSE),
        "must be a character vector of length one"
    )
    expect_error(
        uorfAnalysis(ptn, KozakContext = "medium", plotOut = FALSE),
        "must be either: 'strong'"
    )
    expect_error(
        uorfAnalysis(ptn, onlyUTR5 = "yes", plotOut = FALSE),
        "must be logical"
    )

    testthat::local_mocked_bindings(
        combSeq = function(seqIn) list(tx1 = "TAAAAA", tx2 = "GGGTTT"),
        calc_uORF = function(seqTmp, ext, context, unit) if (grepl("ATG", seqTmp)) 1 else 0,
        .package = "postNet"
    )

    out <- uorfAnalysis(
        ptn,
        startCodon = "ATG",
        KozakContext = "strong",
        onlyUTR5 = FALSE,
        unitOut = "number",
        plotOut = FALSE
    )
    expect_named(out, "uORFs_ATG_strong")
    expect_equal(out$uORFs_ATG_strong, c(gene1 = 1, gene2 = 0))

    testthat::local_mocked_bindings(
        calc_uORF = function(seqTmp, ext, context, unit) if (seqTmp == "AAAATGG") 2 else 1,
        addStats = function(...) NULL,
        roundNice = function(x, direction) 1,
        adjust_ylim = function(lowerLimit, upperLimit) c(0, 1),
        .package = "postNet"
    )
    out_dir <- tempfile("uorf")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)
    uorfAnalysis(
        ptn,
        startCodon = "ATG",
        KozakContext = "strong",
        onlyUTR5 = TRUE,
        unitOut = "number",
        comparisons = list(c(0, 1)),
        plotOut = TRUE,
        pdfName = "uorfplot"
    )
    expect_true(file.exists(file.path(out_dir, "uorfplot_uORFs_strong.pdf")))
})

test_that("goDotplot validates setup and can render pooled and per-list outputs", {
    ptn <- mock_ptn(with_background = TRUE)
    go_slot <- methods::new(
        "postNetGO",
        BP = list(
            comparisonA = methods::new(
                "mockEnrichResult",
                result = data.frame(
                    ID = c("GO:1", "GO:2"),
                    Description = c("term1", "term2"),
                    Count = c(5, 3),
                    Size = c(20, 10),
                    pvalue = c(0.001, 0.002),
                    p.adjust = c(0.01, 0.02),
                    geneID = c("gene1:gene2", "gene2:gene3")
                )
            )
        ),
        CC = NULL,
        MF = NULL,
        KEGG = NULL
    )
    ptn@analysis@GO <- go_slot

    expect_error(goDotplot(ptn, category = "BP", pool = "yes"), "must be logical")
    expect_error(goDotplot(mock_ptn(), category = "BP"), "cannot get a slot")

    out_dir <- tempfile("godot")
    dir.create(out_dir)

    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    expect_null(goDotplot(ptn, category = "BP", pool = TRUE, nCategories = 1, pdfName = "pooled"))
    expect_true(file.exists(file.path(out_dir, "pooled_pooled_GOdotplot_BP_.pdf")))

    expect_null(goDotplot(ptn, category = "BP", pool = FALSE, nCategories = 1, pdfName = "single"))
    expect_true(file.exists(file.path(out_dir, "single_comparisonA_GOdotplot_BP_.pdf")))
})

test_that("goDotplot covers geneRatio, filtering and empty per-list branches", {
    ptn <- mock_ptn(with_background = TRUE)
    go_slot <- methods::new(
        "postNetGO",
        BP = list(
            comparisonA = methods::new(
                "mockEnrichResult",
                result = data.frame(
                    ID = c("GO:1", "GO:2"),
                    Description = c("term1", "term2"),
                    Count = c(5L, 3L),
                    Size = c(20L, 10L),
                    pvalue = c(0.001, 0.002),
                    p.adjust = c(0.01, 0.02),
                    geneID = c("gene1:gene2", "gene2:gene3")
                )
            ),
            comparisonB = methods::new(
                "mockEnrichResult",
                result = data.frame(
                    ID = character(),
                    Description = character(),
                    Count = integer(),
                    Size = integer(),
                    pvalue = numeric(),
                    p.adjust = numeric(),
                    geneID = character()
                )
            )
        ),
        CC = NULL,
        MF = NULL,
        KEGG = NULL
    )
    ptn@analysis@GO <- go_slot

    expect_error(goDotplot(ptn, category = "BP", nCategories = "two"), "nCategories")

    out_dir <- tempfile("godotextra")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    expect_null(
        goDotplot(
            ptn,
            category = "BP",
            pool = TRUE,
            termSel = "GO:2",
            nCategories = 5,
            size = "geneRatio",
            pdfName = "ratio"
        )
    )
    expect_true(file.exists(file.path(out_dir, "ratio_pooled_GOdotplot_BP_.pdf")))

    expect_message(
        goDotplot(
            ptn,
            category = "BP",
            pool = FALSE,
            termSel = "GO:1",
            nCategories = 5,
            size = "geneRatio",
            pdfName = "split"
        ),
        "there are no categories to plot"
    )
    expect_true(file.exists(file.path(out_dir, "split_comparisonA_GOdotplot_BP_.pdf")))
    expect_false(file.exists(file.path(out_dir, "split_comparisonB_GOdotplot_BP_.pdf")))
})

predict.mock_rf <- function(object, newdata, type = c("response", "prob"), ...) {
    type <- match.arg(type)
    if (identical(type, "prob")) {
        out <- cbind(A = rep(0.2, nrow(newdata)), B = rep(0.8, nrow(newdata)))
        return(out)
    }
    factor(rep("B", nrow(newdata)), levels = c("A", "B"))
}

plot.mock_boruta <- function(x, ...) NULL
plot.mock_perf <- function(x, ...) NULL

test_that("featureIntegration covers mocked lm and rf workflows", {
    ptn <- mock_ptn(with_background = TRUE)
    features <- list(
        featureA = c(gene1 = 1, gene2 = 2),
        featureB = c(gene1 = 2, gene2 = 1)
    )

    expect_error(
        featureIntegration(ptn, features = features, analysis_type = "rf", comparisons = NULL),
        "Please provide the desired comparisons"
    )

    testthat::local_mocked_bindings(
        runLM = function(...) {
            methods::new(
                "postNetFeatureIntegration_lm",
                univariateModel = methods::new(
                    "postNetUnivariate",
                    pvalue = c(featureA = 0.01),
                    fdr = c(featureA = 0.02),
                    varianceExplained = c(featureA = 12)
                ),
                stepwiseModel = methods::new(
                    "postNetStepWise",
                    models = list(step1 = "fit"),
                    table = matrix(1, nrow = 1)
                ),
                finalModel = methods::new(
                    "postNetFinalModel",
                    totalVarianceExplained = 20,
                    finalModel = FALSE,
                    table = data.frame(feature = "featureA")
                ),
                selectedFeatures = c(featureA = 1),
                networkGraph = list(nodes = "featureA")
            )
        },
        plotScatterInd = function(...) NULL,
        .package = "postNet"
    )

    lm_out <- featureIntegration(
        ptn,
        features = features,
        analysis_type = "lm",
        comparisons = list(c(0, 1)),
        regOnly = TRUE,
        pdfName = tempfile("lm")
    )
    expect_s4_class(lm_out@analysis@featureIntegration$lm[[1]], "postNetFeatureIntegration_lm")

    testthat::local_mocked_bindings(
        randomForest = function(...) structure(list(importance = matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, dimnames = list(c("featureA", "featureB"), NULL))), class = "mock_rf"),
        importance = function(x) x$importance,
        .package = "randomForest"
    )
    testthat::local_mocked_bindings(
        Boruta = function(...) {
            structure(
                list(
                    x = 1:2,
                    y = 1:2,
                    ImpHistory = matrix(
                        c(1, 2, 3, 4),
                        ncol = 2,
                        dimnames = list(NULL, c("featureA", "shadowMean"))
                    )
                ),
                class = "mock_boruta"
            )
        },
        attStats = function(x) data.frame(
            meanImp = c(5, 1),
            medianImp = c(5, 1),
            minImp = c(4, 0),
            maxImp = c(6, 2),
            normHits = c(1, 0),
            decision = c("Confirmed", "Rejected"),
            row.names = c("featureA", "featureB")
        ),
        .package = "Boruta"
    )
    testthat::local_mocked_bindings(
        prediction = function(predictions, labels) list(predictions = predictions, labels = labels),
        performance = function(pred, measure, x.measure = NULL) {
            if (identical(measure, "auc")) {
                methods::new("mockAuc", y.values = list(0.88))
            } else {
                structure(list(x = c(0, 1), y = c(0, 1)), class = "mock_perf")
            }
        },
        .package = "ROCR"
    )
    testthat::local_mocked_bindings(
        predict = function(object, newdata, type = NULL, ...) {
            if (identical(type, "prob")) {
                return(cbind(A = rep(0.2, nrow(newdata)), B = rep(0.8, nrow(newdata))))
            }
            factor(rep("B", nrow(newdata)), levels = c("A", "B"))
        },
        .package = "stats"
    )
    testthat::local_mocked_bindings(
        confusionMatrix = function(data, reference) list(NULL, NULL, NULL, c(Sensitivity = 0.8, Specificity = 0.9)),
        .package = "caret"
    )
    testthat::local_mocked_bindings(
        pdf = function(...) NULL,
        dev.off = function(...) NULL,
        .package = "grDevices"
    )
    testthat::local_mocked_bindings(
        plotScatterInd = function(...) NULL,
        roundNice = function(x, direction) 10,
        .package = "postNet"
    )

    set.seed(1)
    rf_out <- featureIntegration(
        ptn,
        features = features,
        analysis_type = "rf",
        comparisons = list(c(0, 1)),
        pdfName = tempfile("rf")
    )
    expect_s4_class(rf_out@analysis@featureIntegration$rf[[1]], "postNetFeatureIntegration_rf")
    expect_true("featureA" %in% names(rf_out@analysis@featureIntegration$rf[[1]]@selectedFeatures))
})

test_that("signaturesHeatmap validates inputs and writes a heatmap pdf", {
    ptn <- mock_ptn(with_background = TRUE)
    ptn@dataIn@effect <- c(gene1 = 2, gene2 = -2, gene3 = 0.5, gene4 = -0.5)
    out_dir <- tempfile("sheat")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    expect_error(signaturesHeatmap(ptn, signatureList = list(sigA = "gene1"), unit = "bad"), "can only be 'FDR'")

    testthat::local_mocked_bindings(
        heatmap.2 = function(...) NULL,
        .package = "gplots"
    )

    signaturesHeatmap(
        ptn,
        signatureList = list(sigA = c("gene1", "gene3"), sigB = c("gene2", "gene4")),
        unit = "FDR",
        pdfName = "sigheat"
    )
    expect_true(file.exists(file.path(out_dir, "sigheat_heatmap.pdf")))
})

test_that("plotFeaturesMap validates inputs and stores mocked UMAP coordinates", {
    ptn <- mock_ptn(with_background = TRUE)

    expect_error(plotFeaturesMap(ptn, featSel = "missing"), "Please provide a character vector of features of interest")
    ptn@features <- data.frame(featureA = c(1, 2), featureB = c(2, 1), row.names = c("gene1", "gene2"))
    expect_error(plotFeaturesMap(ptn, featSel = c("featureA", "featureB"), regOnly = TRUE, comparisons = list(c(0, 1), c(0, 1))), "please run each comparison separately")
    expect_error(plotFeaturesMap(ptn, featSel = c("featureA", "featureB"), regOnly = FALSE, remExtreme = 2), "between 0 and 1")

    testthat::local_mocked_bindings(
        umap = function(data, n_components = 2) {
            list(layout = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE, dimnames = list(rownames(data), NULL)))
        },
        .package = "umap"
    )
    testthat::local_mocked_bindings(
        plot_fmap = function(fMap, colVec, remExtreme = NULL, name) {
            list(mainPlot = name, legend = list(name))
        },
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        grid.arrange = function(...) NULL,
        .package = "gridExtra"
    )

    out_dir <- tempfile("fmap")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    out <- plotFeaturesMap(
        ptn,
        regOnly = FALSE,
        featSel = c("featureA", "featureB"),
        featCol = "featureA",
        pdfName = "fmap"
    )
    expect_true(file.exists(file.path(out_dir, "fmap_featureA_featureUMAP.pdf")))
    expect_equal(colnames(out@analysis@featureIntegration$featuresMap), c("UMAP1", "UMAP2"))

    ptn@features <- data.frame(
        featureA = c(1, 2, 3),
        featureB = c(2, 1, 3),
        binary = c(0, 1, 0),
        row.names = c("gene1", "gene2", "gene3")
    )
    testthat::local_mocked_bindings(
        umap = function(data, n_components = 2) {
            expect_true("binary" %in% colnames(data))
            list(layout = matrix(c(0, 0, 1, 1, 2, 2), ncol = 2, byrow = TRUE, dimnames = list(rownames(data), NULL)))
        },
        .package = "umap"
    )
    out_scaled <- plotFeaturesMap(
        ptn,
        regOnly = FALSE,
        featSel = c("featureA", "featureB", "binary"),
        featCol = c("featureA", "binary"),
        remBinary = FALSE,
        scaled = TRUE,
        centered = TRUE,
        remExtreme = 0.1,
        pdfName = "fmapscaled"
    )
    expect_true(file.exists(file.path(out_dir, "fmapscaled_featureA_featureUMAP.pdf")))
    expect_true(file.exists(file.path(out_dir, "fmapscaled_binary_featureUMAP.pdf")))
    expect_equal(colnames(out_scaled@analysis@featureIntegration$featuresMap), c("UMAP1", "UMAP2"))
})

test_that("plotFeaturesMap covers additional validation and regulated-only branch", {
    ptn <- mock_ptn(with_background = TRUE)
    ptn@features <- data.frame(
        featureA = c(1, 2, 3),
        featureB = c(2, 1, 3),
        binary = c(0, 1, 0),
        row.names = c("gene1", "gene2", "gene3")
    )

    expect_error(plotFeaturesMap(ptn, featSel = c("featureA", "featureB"), featCol = "missing"), "that will be overlaid")
    expect_error(plotFeaturesMap(ptn, featSel = c("featureA", "featureB"), regOnly = "yes"), "regOnly")
    expect_error(plotFeaturesMap(ptn, featSel = c("featureA", "featureB"), regOnly = FALSE, scaled = "yes"), "scaled")
    expect_error(plotFeaturesMap(ptn, featSel = c("featureA", "featureB"), regOnly = FALSE, centered = "yes"), "centered")
    expect_error(plotFeaturesMap(ptn, featSel = c("featureA", "featureB"), regOnly = FALSE, remBinary = "yes"), "remBinary")
    expect_error(plotFeaturesMap(ptn, featSel = c("featureA", "featureB"), regOnly = TRUE, comparisons = "bad"), "must be a list of numeric")

    ptn_no_bg <- ptn
    ptn_no_bg@dataIn@background <- NULL
    expect_error(
        plotFeaturesMap(ptn_no_bg, featSel = c("featureA", "featureB"), regOnly = TRUE, comparisons = list(c(0, 1))),
        "0 always denotes the background"
    )

    testthat::local_mocked_bindings(
        umap = function(data, n_components = 2) {
            expect_false("binary" %in% colnames(as.data.frame(data)))
            list(layout = matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE, dimnames = list(rownames(data), NULL)))
        },
        .package = "umap"
    )
    testthat::local_mocked_bindings(
        plot_fmap = function(fMap, colVec, remExtreme = NULL, name) {
            list(mainPlot = name, legend = list(name))
        },
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        grid.arrange = function(...) NULL,
        .package = "gridExtra"
    )

    out_dir <- tempfile("fmapreg")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    out_reg <- plotFeaturesMap(
        ptn,
        regOnly = TRUE,
        comparisons = list(c(0, 1)),
        featSel = c("featureA", "featureB", "binary"),
        featCol = c("featureA", "featureB"),
        remBinary = TRUE,
        scaled = FALSE,
        pdfName = "fmapreg"
    )
    expect_equal(colnames(out_reg@analysis@featureIntegration$featuresMap), c("UMAP1", "UMAP2"))
    expect_true(file.exists(file.path(out_dir, "fmapreg_featureA_featureUMAP.pdf")))
    expect_true(file.exists(file.path(out_dir, "fmapreg_featureB_featureUMAP.pdf")))
})

test_that("plotSignatures validates inputs and writes a signature pdf", {
    ptn <- mock_ptn(with_background = TRUE)
    ptn@dataIn@effect <- c(gene1 = 2, gene2 = -2, gene3 = 0.5, gene4 = -0.5)
    sigs <- list(sigA = "gene1", sigB = "gene2")

    expect_error(plotSignatures(ptn, signatureList = sigs, dataName = NULL, generalName = "test"), "Please provide an input for 'generalName'")
    expect_error(plotSignatures(ptn, signatureList = sigs, dataName = "data", generalName = "test", tableCex = "big"), "numeric value for 'tableCex'")

    out_dir <- tempfile("psig")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    plotSignatures(
        ptn,
        signatureList = sigs,
        dataName = "data",
        generalName = "test",
        signature_colours = c("red", "blue"),
        pdfName = "sigplot"
    )
    expect_true(file.exists(file.path(out_dir, "sigplot_data_data_signature_test.pdf")))

    testthat::local_mocked_bindings(
        addtable2plot = function(...) NULL,
        .package = "plotrix"
    )
    expect_message(
        plotSignatures(
            ptn,
            signatureList = list(sigA = c("gene1", "gene2"), sigB = c("gene2", "gene3")),
            dataName = "data2",
            generalName = "overlap",
            signature_colours = NULL,
            xlim = c(-3, 3),
            pdfName = "sigplot2"
        ),
        "overlap between signatures"
    )
    expect_true(file.exists(file.path(out_dir, "sigplot2_data_data2_signature_overlap.pdf")))
})

test_that("gseaPlot validates inputs and writes a term plot", {
    ptn <- mock_ptn(with_background = TRUE, with_gsea = TRUE)

    expect_error(gseaPlot(ptn, termNames = "setA", gseaParam = "x"), "numeric values for 'gseaParam'")

    testthat::local_mocked_bindings(
        calcGseaStat = function(stats, selectedStats, returnAllExtremes = TRUE) {
            list(bottoms = c(0, -0.2), tops = c(0.3, 0.5))
        },
        .package = "fgsea"
    )

    out_dir <- tempfile("gseaplot")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    gseaPlot(ptn, termNames = "setA", pdfName = "gplot")
    expect_true(file.exists(file.path(out_dir, "gplot_gsea_setA.pdf")))
})

test_that("codonUsage validates inputs and stores mocked codon analysis", {
    ptn <- mock_ptn(with_background = TRUE)
    methods::slot(ptn@annot, "CDS") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("ATGAAAAAA", "ATGAAAGGG")
    )

    expect_error(codonUsage(ptn, analysis = "AA", comparisons = NULL), "pairs of comparisons must be specified")
    expect_error(codonUsage(ptn, analysis = "bad", comparisons = list(c(0, 1))), "options are 'codon' or 'AA'")
    expect_error(codonUsage(ptn, annotType = "ccds", sourceSeq = "bad", analysis = "AA", comparisons = list(c(0, 1))), "sourceSeq")
    expect_error(codonUsage(ptn, analysis = "AA", codonN = "one", comparisons = list(c(0, 1))), "integer value for 'codonN'")
    expect_error(codonUsage(ptn, analysis = "AA", rem5 = "yes", comparisons = list(c(0, 1))), "input for 'rem5' must be logical")
    expect_error(codonUsage(ptn, analysis = "AA", subregion = c(1, 2), comparisons = list(c(0, 1))), "input for 'subregion' must be an integer")
    expect_error(codonUsage(ptn, analysis = "AA", subregion = 2, subregionSel = "keep", comparisons = list(c(0, 1))), "must be either 'select' or 'exclude'")

    testthat::local_mocked_bindings(
        chisq.test = function(x) list(stdres = matrix(c(3, -3), nrow = 1, dimnames = list("A", c("background", "comparisonA"))), p.value = 0.001),
        .package = "stats"
    )
    testthat::local_mocked_bindings(
        codonCount = function(seq, gene, codonN = 1) {
            data.frame(
                geneID = gene,
                codon = c("AAA", "GGG"),
                AA = c("K", "G"),
                count = c(10, 8),
                frequency = c(0.6, 0.4),
                AACountPerGene = c(10, 8),
                relative_frequency = c(1, 1)
            )
        },
        plotPostNet = function(...) NULL,
        statOnDf = function(df, regs, analysis) {
            nm <- rownames(df)[1]
            if (is.null(nm) || is.na(nm) || identical(nm, "")) {
                nm <- as.character(df[[1]][1])
            }
            stats::setNames(4, nm)
        },
        roundNice = function(x, direction) 1,
        .package = "postNet"
    )

    out_dir <- tempfile("codon")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    out <- codonUsage(
        ptn,
        analysis = "AA",
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = FALSE,
        pdfName = "codon"
    )

    expect_s4_class(out@analysis@codons, "postNetCodons")
    expect_equal(length(out@analysis@codons@codonSelection), 0)
    expect_equal(nrow(ptn_codonAnalysis(out)), 4)

    testthat::local_mocked_bindings(
        codonCount = function(seq, gene, codonN = 1) {
            data.frame(
                geneID = gene,
                codon = c("AAA", "GGG"),
                AA = c("K", "G"),
                count = c(12, 6),
                frequency = c(0.6, 0.4),
                AACountPerGene = c(12, 6),
                relative_frequency = c(1, 1)
            )
        },
        statOnDf = function(df, regs, analysis) stats::setNames(c(4, 0.5), c("K", "G")),
        roundNice = function(x, direction) 2,
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        chisq.test = function(x) list(
            stdres = matrix(
                c(4, -4, 0.5, -0.5),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(c("background", "comparisonA"), c("K", "G"))
            ),
            p.value = 0.001
        ),
        .package = "stats"
    )
    testthat::local_mocked_bindings(
        heatmap.2 = function(...) NULL,
        .package = "gplots"
    )
    out_aa_plot <- codonUsage(
        ptn,
        analysis = "AA",
        comparisons = list(c(0, 1)),
        plotHeatmap = TRUE,
        rem5 = FALSE,
        pdfName = "codonAA"
    )
    expect_s4_class(out_aa_plot@analysis@codons, "postNetCodons")
    expect_true(sum(grepl("^codonAA_.*\\.pdf$", list.files(out_dir))) >= 1)

    ptn@species <- "human"
    testthat::local_mocked_bindings(
        codonCount = function(seq, gene, codonN = 1) {
            data.frame(
                geneID = gene,
                codon = c("AAA", "AAG", "GGA"),
                AA = c("K", "K", "G"),
                count = c(10, 8, 6),
                frequency = c(0.4, 0.35, 0.25),
                AACountPerGene = c(18, 18, 6),
                relative_frequency = c(10 / 18, 8 / 18, 1)
            )
        },
        get_reference_data = function(file) data.frame(
            external_gene_name = c("gene1", "gene2", "gene3"),
            CAI = c(0.1, 0.2, 0.3),
            CBI = c(0.1, 0.2, 0.3),
            Fop = c(0.1, 0.2, 0.3),
            tAI = c(0.1, 0.2, 0.3),
            L_aa = c(1, 2, 3)
        ),
        plotPostNet = function(...) NULL,
        statOnDf = function(df, regs, analysis) stats::setNames(c(4, 0.5), c("AAA", "AAG")),
        roundNice = function(x, direction) 1,
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        read.delim = function(...) data.frame(
            external_gene_name = c("gene1", "gene2", "gene3"),
            CAI = c(0.1, 0.2, 0.3),
            CBI = c(0.1, 0.2, 0.3),
            Fop = c(0.1, 0.2, 0.3),
            tAI = c(0.1, 0.2, 0.3),
            L_aa = c(1, 2, 3)
        ),
        .package = "utils"
    )
    testthat::local_mocked_bindings(
        chisq.test = function(x) list(
            stdres = matrix(c(4, -4, 0.5, -0.5), nrow = 1, dimnames = list("A", c("AAA", "AAG", "GGA", "GGG"))),
            p.value = 0.001
        ),
        .package = "stats"
    )

    out_codon <- codonUsage(
        ptn,
        analysis = "codon",
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = FALSE,
        pdfName = "codon2"
    )
    expect_s4_class(out_codon@analysis@codons, "postNetCodons")

    testthat::local_mocked_bindings(
        codonCount = function(seq, gene, codonN = 1) {
            data.frame(
                geneID = gene,
                codon = c("AAA", "AAG", "GGA", "GGG"),
                AA = c("K", "K", "G", "G"),
                count = c(12, 8, 7, 5),
                frequency = c(0.4, 0.3, 0.2, 0.1),
                AACountPerGene = c(20, 20, 12, 12),
                relative_frequency = c(0.6, 0.4, 7 / 12, 5 / 12)
            )
        },
        statOnDf = function(df, regs, analysis) {
            stats::setNames(c(4, 2, 0.5, 0.25), c("AAA", "AAG", "GGA", "GGG"))
        },
        roundNice = function(x, direction) 2,
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        chisq.test = function(x) list(
            stdres = matrix(
                c(4, -4, 3.5, -3.5, 0.5, -0.5, 0.25, -0.25),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(c("background", "comparisonA"), c("AAA", "AAG", "GGA", "GGG"))
            ),
            p.value = 0.001
        ),
        .package = "stats"
    )
    testthat::local_mocked_bindings(
        read.delim = function(...) data.frame(
            external_gene_name = c("gene1", "gene2", "gene3"),
            CAI = c(0.1, 0.2, 0.3),
            CBI = c(0.1, 0.2, 0.3),
            Fop = c(0.1, 0.2, 0.3),
            tAI = c(0.1, 0.2, 0.3),
            L_aa = c(1, 2, 3)
        ),
        .package = "utils"
    )
    testthat::local_mocked_bindings(
        heatmap.2 = function(...) NULL,
        .package = "gplots"
    )
    out_plot <- codonUsage(
        ptn,
        analysis = "codon",
        comparisons = list(c(0, 1)),
        plotHeatmap = TRUE,
        rem5 = FALSE,
        subregion = 6,
        subregionSel = "select",
        pdfName = "codon3"
    )
    expect_s4_class(out_plot@analysis@codons, "postNetCodons")
    expect_true(length(list.files(out_dir, pattern = "^codon3_.*pdf$")) >= 5)
})

test_that("codonUsage supports CCDS annotations from loaded and created sources", {
    ptn <- mock_ptn(with_background = TRUE)
    ptn@species <- "human"

    ptn_bad_species <- ptn
    ptn_bad_species@species <- "rat"
    expect_error(
        codonUsage(
            ptn_bad_species,
            annotType = "ccds",
            sourceSeq = "load",
            analysis = "AA",
            comparisons = list(c(0, 1)),
            plotHeatmap = FALSE,
            rem5 = FALSE
        ),
        "Please specify a species"
    )

    testthat::local_mocked_bindings(
        get_reference_data = function(file) data.frame(
            id = c("CCDS1", "CCDS2"),
            geneID = c("gene1", "gene2"),
            CDS_seq = c("ATGAAATTT", "ATGGGGCCC")
        ),
        codonCount = function(seq, gene, codonN = 1) {
            data.frame(
                geneID = gene,
                codon = c("AAA", "GGG"),
                AA = c("K", "G"),
                count = c(6, 4),
                frequency = c(0.6, 0.4),
                AACountPerGene = c(6, 4),
                relative_frequency = c(1, 1)
            )
        },
        statOnDf = function(df, regs, analysis) stats::setNames(c(4, 0.5), c("AAA", "GGG")),
        roundNice = function(x, direction) 1,
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        chisq.test = function(x) list(
            stdres = matrix(c(4, -4), nrow = 1, dimnames = list("background", c("AAA", "GGG"))),
            p.value = 0.001
        ),
        .package = "stats"
    )
    out_load <- codonUsage(
        ptn,
        annotType = "ccds",
        sourceSeq = "load",
        analysis = "AA",
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = FALSE
    )
    expect_s4_class(out_load@annot@CCDS, "postNetRegion")
    expect_equal(ptn_geneID(out_load, region = "CCDS"), c("gene1", "gene2"))

    testthat::local_mocked_bindings(
        download.file = function(url, destfile, ...) file.create(destfile),
        read.delim = function(file, ...) data.frame(
            ccds_id = c("CCDS1.1", "CCDS2.1"),
            gene = c("gene1", "gene2")
        ),
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        download.file = function(url, destfile, ...) file.create(destfile),
        .package = "utils"
    )
    testthat::local_mocked_bindings(
        unlink = function(x, ...) TRUE,
        .package = "base"
    )
    testthat::local_mocked_bindings(
        read.delim = function(file, ...) data.frame(
            ccds_id = c("CCDS1.1", "CCDS2.1"),
            gene = c("gene1", "gene2")
        ),
        .package = "utils"
    )
    testthat::local_mocked_bindings(
        gunzip = function(filename, remove = TRUE) sub("\\.gz$", "", filename),
        .package = "R.utils"
    )
    testthat::local_mocked_bindings(
        read.fasta = function(file, seqtype = "AA") {
            list(
                "CCDS1.1" = strsplit("ATGAAATTT", "")[[1]],
                "CCDS2.1" = strsplit("ATGGGGCCC", "")[[1]]
            )
        },
        .package = "seqinr"
    )
    out_create <- codonUsage(
        ptn,
        annotType = "ccds",
        sourceSeq = "create",
        analysis = "AA",
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = FALSE
    )
    expect_s4_class(out_create@annot@CCDS, "postNetRegion")
    expect_equal(length(ptn_geneID(out_create, region = "CCDS")), 2)
})

test_that("codonUsage reports contingency-table edge cases clearly", {
    ptn <- mock_ptn(with_background = TRUE)
    methods::slot(ptn@annot, "CDS") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("ATGAAAAAA", "ATGAAAGGG")
    )

    testthat::local_mocked_bindings(
        codonCount = function(seq, gene, codonN = 1) {
            data.frame(
                geneID = gene,
                codon = c("AAA", "GGG"),
                AA = c("K", "G"),
                count = c(0, 0),
                frequency = c(0, 0),
                AACountPerGene = c(0, 0),
                relative_frequency = c(0, 0)
            )
        },
        .package = "postNet"
    )
    expect_error(
        codonUsage(
            ptn,
            analysis = "AA",
            comparisons = list(c(0, 1)),
            plotHeatmap = FALSE,
            rem5 = FALSE
        ),
        "must have at least one positive entry"
    )

})

test_that("codonUsage covers low-count warnings", {
    ptn <- mock_ptn(with_background = TRUE)
    methods::slot(ptn@annot, "CDS") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("ATGAAAAAA", "ATGAAAGGG")
    )

    testthat::local_mocked_bindings(
        codonCount = function(seq, gene, codonN = 1) {
            data.frame(
                geneID = gene,
                codon = c("AAAAAA", "AAAGGG"),
                AA = c("KK", "KG"),
                count = c(2, 1),
                frequency = c(0.67, 0.33),
                AACountPerGene = c(2, 1),
                relative_frequency = c(1, 1)
            )
        },
        statOnDf = function(df, regs, analysis) stats::setNames(c(4, 2), c("AAAAAA", "AAAGGG")),
        .package = "postNet"
    )

    expect_warning(
        codonUsage(
            ptn,
            analysis = "codon",
            codonN = 2,
            comparisons = list(c(0, 1)),
            plotHeatmap = FALSE,
            rem5 = FALSE
        ),
        "some counts are lower than 5"
    )
})

test_that("codonUsage covers codonN > 1 selection branch with significant di-codons", {
    ptn <- mock_ptn(with_background = TRUE)
    methods::slot(ptn@annot, "CDS") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("ATGAAAAAA", "ATGAAAGGG")
    )

    testthat::local_mocked_bindings(
        codonCount = function(seq, gene, codonN = 1) {
            if (gene == "gene1") {
                data.frame(
                    geneID = gene,
                    codon = c("AAAAAA", "AAAGGG"),
                    AA = c("KK", "KG"),
                    count = c(30, 8),
                    frequency = c(0.79, 0.21),
                    AACountPerGene = c(30, 8),
                    relative_frequency = c(1, 1)
                )
            } else {
                data.frame(
                    geneID = gene,
                    codon = c("AAAAAA", "AAAGGG"),
                    AA = c("KK", "KG"),
                    count = c(6, 24),
                    frequency = c(0.2, 0.8),
                    AACountPerGene = c(6, 24),
                    relative_frequency = c(1, 1)
                )
            }
        },
        .package = "postNet"
    )

    out_dicodon <- codonUsage(
        ptn,
        analysis = "codon",
        codonN = 2,
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = TRUE
    )
    expect_s4_class(out_dicodon@analysis@codons, "postNetCodons")
    expect_true(is.list(out_dicodon@analysis@codons@codonSelection))
})

test_that("codonUsage covers additional mouse and message branches", {
    ptn <- mock_ptn(with_background = TRUE)
    methods::slot(ptn@annot, "CDS") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("ATGAAAAAA", "ATGAAAGGG")
    )

    expect_error(
        codonUsage(ptn, annotType = "bad", analysis = "AA", comparisons = list(c(0, 1))),
        "source of annotation"
    )

    ptn_no_bg <- mock_ptn(with_background = FALSE)
    methods::slot(ptn_no_bg@annot, "CDS") <- methods::slot(ptn@annot, "CDS")
    expect_error(
        codonUsage(ptn_no_bg, analysis = "AA", comparisons = list(c(0, 1))),
        "0 always denotes the background"
    )
    expect_error(
        codonUsage(ptn, analysis = "AA", comparisons = "bad"),
        "must be a list of numeric"
    )

    testthat::local_mocked_bindings(
        codonCount = function(seq, gene, codonN = 1) {
            data.frame(
                geneID = gene,
                codon = c("AAA", "GGG"),
                AA = c("K", "G"),
                count = c(10, 8),
                frequency = c(0.6, 0.4),
                AACountPerGene = c(10, 8),
                relative_frequency = c(1, 1)
            )
        },
        .package = "postNet"
    )

    ptn_bad_species <- ptn
    ptn_bad_species@species <- "rat"
    expect_error(
        codonUsage(
            ptn_bad_species,
            analysis = "codon",
            comparisons = list(c(0, 1)),
            plotHeatmap = FALSE,
            rem5 = FALSE
        ),
        "codon usage indexes"
    )

    testthat::local_mocked_bindings(
        chisq.test = function(x) list(
            stdres = matrix(
                c(0.1, -0.1, 0.05, -0.05),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(c("background", "comparisonA"), c("K", "G"))
            ),
            p.value = 0.001
        ),
        .package = "stats"
    )
    expect_message(
        no_signif <- codonUsage(
            ptn,
            analysis = "AA",
            comparisons = list(c(0, 1)),
            plotHeatmap = FALSE,
            rem5 = FALSE
        ),
        "There are no significant codons for the comparison 1"
    )
    expect_s4_class(no_signif@analysis@codons, "postNetCodons")

    ptn_mouse <- ptn
    ptn_mouse@species <- "mouse"
    testthat::local_mocked_bindings(
        get_reference_data = function(file) {
            if (identical(file, "mouseDB_ccds.txt.gz")) {
                data.frame(
                    id = c("CCDSm1", "CCDSm2"),
                    geneID = c("gene1", "gene2"),
                    CDS_seq = c("ATGAAATTT", "ATGGGGCCC")
                )
            } else {
                data.frame(
                    external_gene_name = c("gene1", "gene2", "gene3"),
                    CAI = c(0.1, 0.2, 0.3),
                    CBI = c(0.1, 0.2, 0.3),
                    Fop = c(0.1, 0.2, 0.3),
                    tAI = c(0.1, 0.2, 0.3),
                    L_aa = c(1, 2, 3)
                )
            }
        },
        statOnDf = function(df, regs, analysis) stats::setNames(c(4, 0.5), c("AAA", "GGG")),
        roundNice = function(x, direction) 1,
        plotPostNet = function(...) NULL,
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        read.delim = function(...) data.frame(
            external_gene_name = c("gene1", "gene2", "gene3"),
            CAI = c(0.1, 0.2, 0.3),
            CBI = c(0.1, 0.2, 0.3),
            Fop = c(0.1, 0.2, 0.3),
            tAI = c(0.1, 0.2, 0.3),
            L_aa = c(1, 2, 3)
        ),
        download.file = function(url, destfile, ...) file.create(destfile),
        .package = "utils"
    )
    testthat::local_mocked_bindings(
        read.delim = function(file, ...) data.frame(
            ccds_id = c("CCDSm1.1", "CCDSm2.1"),
            gene = c("gene1", "gene2")
        ),
        download.file = function(url, destfile, ...) file.create(destfile),
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        unlink = function(x, ...) TRUE,
        .package = "base"
    )
    testthat::local_mocked_bindings(
        gunzip = function(filename, remove = TRUE) sub("\\.gz$", "", filename),
        .package = "R.utils"
    )
    testthat::local_mocked_bindings(
        read.fasta = function(file, seqtype = "AA") {
            list(
                "CCDSm1.1" = strsplit("ATGAAATTT", "")[[1]],
                "CCDSm2.1" = strsplit("ATGGGGCCC", "")[[1]]
            )
        },
        .package = "seqinr"
    )
    testthat::local_mocked_bindings(
        chisq.test = function(x) list(
            stdres = matrix(
                c(4, -4, 0.5, -0.5),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(c("background", "comparisonA"), c("AAA", "GGG"))
            ),
            p.value = 0.001
        ),
        .package = "stats"
    )

    out_mouse_load <- codonUsage(
        ptn_mouse,
        annotType = "ccds",
        sourceSeq = "load",
        analysis = "AA",
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = FALSE
    )
    expect_s4_class(out_mouse_load@annot@CCDS, "postNetRegion")

    out_mouse_create <- codonUsage(
        ptn_mouse,
        annotType = "ccds",
        sourceSeq = "create",
        analysis = "AA",
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = FALSE
    )
    expect_s4_class(out_mouse_create@annot@CCDS, "postNetRegion")

    out_mouse_codon <- codonUsage(
        ptn_mouse,
        analysis = "codon",
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = FALSE,
        pdfName = "mouseCodon"
    )
    expect_s4_class(out_mouse_codon@analysis@codons, "postNetCodons")
})

test_that("codonUsage covers codon and AA label-threshold plotting branches", {
    ptn <- mock_ptn(with_background = TRUE)
    methods::slot(ptn@annot, "CDS") <- mock_region(
        ids = c("tx1", "tx2"),
        gene_ids = c("gene1", "gene2"),
        sequences = c("ATGAAAAAA", "ATGAAAGGG")
    )

    testthat::local_mocked_bindings(
        codonCount = function(seq, gene, codonN = 1) {
            data.frame(
                geneID = gene,
                codon = c("AAA", "AAG", "GGA", "GGG"),
                AA = c("K", "K", "G", "G"),
                count = c(20, 16, 12, 8),
                frequency = c(0.4, 0.3, 0.2, 0.1),
                AACountPerGene = c(36, 36, 20, 20),
                relative_frequency = c(20 / 36, 16 / 36, 12 / 20, 8 / 20)
            )
        },
        get_reference_data = function(file) data.frame(
            external_gene_name = c("gene1", "gene2", "gene3"),
            CAI = c(0.1, 0.2, 0.3),
            CBI = c(0.1, 0.2, 0.3),
            Fop = c(0.1, 0.2, 0.3),
            tAI = c(0.1, 0.2, 0.3),
            L_aa = c(1, 2, 3)
        ),
        plotPostNet = function(...) NULL,
        roundNice = function(x, direction) 2,
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        read.delim = function(...) data.frame(
            external_gene_name = c("gene1", "gene2", "gene3"),
            CAI = c(0.1, 0.2, 0.3),
            CBI = c(0.1, 0.2, 0.3),
            Fop = c(0.1, 0.2, 0.3),
            tAI = c(0.1, 0.2, 0.3),
            L_aa = c(1, 2, 3)
        ),
        .package = "utils"
    )
    testthat::local_mocked_bindings(
        pdf = function(...) NULL,
        dev.off = function(...) NULL,
        .package = "grDevices"
    )

    ptn@species <- "human"
    testthat::local_mocked_bindings(
        chisq.test = function(x) list(
            stdres = matrix(
                c(4, 4, 4, 4, -4, -4, -4, -4),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(c("background", "comparisonA"), c("AAA", "AAG", "GGA", "GGG"))
            ),
            p.value = 0.001
        ),
        .package = "stats"
    )
    testthat::local_mocked_bindings(
        statOnDf = function(df, regs, analysis) {
            stats::setNames(c(0.25, 0.125, 0.0625, 0.03125), c("AAA", "AAG", "GGA", "GGG"))
        },
        .package = "postNet"
    )
    out_codon_msg <- codonUsage(
        ptn,
        analysis = "codon",
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = FALSE
    )
    expect_s4_class(out_codon_msg@analysis@codons, "postNetCodons")

    testthat::local_mocked_bindings(
        chisq.test = function(x) list(
            stdres = matrix(
                c(4, 4, -4, -4),
                nrow = 2,
                byrow = TRUE,
                dimnames = list(c("background", "comparisonA"), c("K", "G"))
            ),
            p.value = 0.001
        ),
        .package = "stats"
    )
    testthat::local_mocked_bindings(
        statOnDf = function(df, regs, analysis) stats::setNames(c(0.25, 0.125), c("K", "G")),
        .package = "postNet"
    )
    out_aa_msg <- codonUsage(
        ptn,
        analysis = "AA",
        comparisons = list(c(0, 1)),
        plotHeatmap = FALSE,
        rem5 = FALSE
    )
    expect_s4_class(out_aa_msg@analysis@codons, "postNetCodons")
})

test_that("rfPred validates inputs and runs with a mocked rf model", {
    ptn <- mock_ptn(with_feature_models = TRUE)
    ptn@analysis@featureIntegration$rf[[1]]@finalModel <- structure(list(), class = "mock_rf")
    ptn@analysis@featureIntegration$rf[[1]]@selectedFeatures <- c(featureA = 1, featureB = 1)

    pred_features <- list(
        featureA = c(gene1 = 1, gene2 = 2),
        featureB = c(gene1 = 2, gene2 = 1)
    )
    pred_gene_list <- list(groupA = "gene1", groupB = "gene2")

    expect_error(rfPred(mock_ptn(), comparison = 1, predGeneList = pred_gene_list, predFeatures = pred_features), "Random Forest implementation first")
    expect_error(
        rfPred(
            ptn,
            comparison = 1,
            predGeneList = pred_gene_list,
            predFeatures = list(
                featureB = c(gene1 = 1, gene2 = 2),
                featureC = c(gene1 = 2, gene2 = 1)
            )
        ),
        "are missing in the 'predFeatures' object"
    )

    testthat::local_mocked_bindings(
        predict = function(object, newdata, type = NULL, ...) {
            if (identical(type, "prob")) {
                return(cbind(A = rep(0.2, nrow(newdata)), B = rep(0.8, nrow(newdata))))
            }
            factor(rep("B", nrow(newdata)), levels = c("A", "B"))
        },
        .package = "stats"
    )
    testthat::local_mocked_bindings(
        prediction = function(predictions, labels) list(predictions = predictions, labels = labels),
        performance = function(pred, measure, x.measure = NULL) {
            structure(list(x = c(0, 1), y = c(0, 1)), class = "mock_perf")
        },
        .package = "ROCR"
    )
    testthat::local_mocked_bindings(
        confusionMatrix = function(data, reference) list(NULL, NULL, c(Accuracy = 0.85), c(Sensitivity = 0.8, Specificity = 0.9)),
        .package = "caret"
    )

    out_dir <- tempfile("rfpred")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    out <- rfPred(
        ptn,
        comparison = 1,
        predGeneList = pred_gene_list,
        predFeatures = pred_features,
        pdfName = "rfpred"
    )
    expect_s4_class(out, "postNetData")
    expect_true(file.exists(file.path(out_dir, "rfpred_rocr.pdf")))
})

test_that("plotSignatures_ads validates inputs and writes a pdf", {
    ads_min <- mock_ads()
    sigs <- list(sigA = "gene1", sigB = "gene2")

    expect_error(plotSignatures_ads(ads_min, contrast = 3, dataName = "data", signatureList = sigs, generalName = "sig", signature_colours = c("red", "blue")), "corresponding to the number of the anota2seq contrast")
    expect_error(plotSignatures_ads(ads_min, contrast = 1, dataName = "data", effects_names = c("a", "b", "c"), signatureList = sigs, generalName = "sig", signature_colours = c("red", "blue")), "There must be 4 effect names provided")
    expect_error(plotSignatures_ads(ads_min, contrast = 1, dataName = "data", signatureList = sigs, generalName = NULL, signature_colours = c("red", "blue")), "Please provide an input for 'generalName'")
    expect_error(plotSignatures_ads(ads_min, contrast = 1, dataName = "data", signatureList = sigs, generalName = "sig", signature_colours = "red"), "signature_colours")
    expect_error(plotSignatures_ads(ads_min, contrast = 1, dataName = "data", signatureList = sigs, generalName = "sig", signature_colours = c("red", "blue"), tableCex = "large"), "tableCex")
    expect_error(plotSignatures_ads(ads_min, contrast = 1, dataName = "data", signatureList = sigs, generalName = "sig", signature_colours = c("red", "blue"), scatterXY = "wide"), "scatterXY")

    testthat::skip_if_not("totalmRNA" %in% names(methods::getSlots("Anota2seqDataSet")))
    ads <- mock_ads_plot()

    expect_error(plotSignatures_ads(ads, contrast = 2, dataName = "data", signatureList = sigs, generalName = "sig", signature_colours = c("red", "blue")), "corresponding to the number of the anota2seq contrast")
    expect_error(plotSignatures_ads(ads, contrast = 1, dataName = NULL, signatureList = sigs, generalName = "sig", signature_colours = c("red", "blue")), "Please provide an input for 'dataName'")

    out_dir <- tempfile("psa")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    plotSignatures_ads(
        ads,
        contrast = 1,
        dataName = "data",
        signatureList = sigs,
        generalName = "sig",
        signature_colours = c("red", "blue"),
        pdfName = "adsplot"
    )
    expect_true(file.exists(file.path(out_dir, "adsplot_data_data_signature_sig.pdf")))

    testthat::local_mocked_bindings(
        addtable2plot = function(...) NULL,
        .package = "plotrix"
    )
    expect_message(
        plotSignatures_ads(
            ads,
            contrast = 1,
            dataName = "data",
            signatureList = list(sigA = c("gene1", "gene2"), sigB = c("gene2", "gene3")),
            generalName = "sigOverlap",
            signature_colours = c("red", "blue"),
            xlim = NULL,
            scatterXY = 3,
            pdfName = "adsplot2"
        ),
        "overlap between signatures"
    )
    expect_true(file.exists(file.path(out_dir, "adsplot2_data_data_signature_sigOverlap.pdf")))
})

test_that("postNetStart can build postNetData from mocked annotation sources", {
    gene_list <- list(groupA = c("gene1", "gene2", "gene3", "gene4"))
    effect_measure <- c(gene1 = 1, gene2 = 2, gene3 = -1, gene4 = -2)

    custom_file <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            id = c("tx1", "tx2", "tx3", "tx4"),
            geneID = c("gene1", "gene2", "gene3", "gene4"),
            UTR5_seq = c("AAA", "CCC", "GGG", "TTT"),
            CDS_seq = c("ATGAAA", "ATGCCC", "ATGGGG", "ATGTTT"),
            UTR3_seq = c("TTT", "AAA", "CCC", "GGG")
        ),
        file = custom_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    ptn_custom <- postNetStart(
        geneList = gene_list,
        geneListcolours = "red",
        effectMeasure = effect_measure,
        source = "custom",
        customFile = custom_file,
        selection = "random"
    )
    expect_s4_class(ptn_custom, "postNetData")
    expect_equal(ptn_geneList(ptn_custom), gene_list)

    testthat::local_mocked_bindings(
        get_reference_data = function(file) data.frame(
            id = c("tx1", "tx2", "tx3", "tx4"),
            geneID = c("gene1", "gene2", "gene3", "gene4"),
            UTR5_seq = c("AAA", "CCC", "GGG", "TTT"),
            CDS_seq = c("ATGAAA", "ATGCCC", "ATGGGG", "ATGTTT"),
            UTR3_seq = c("TTT", "AAA", "CCC", "GGG")
        ),
        anota2seqGetDirectedRegulations = function(ads) {
            list(list(
                translationUp = c("gene1", "gene2", "gene3", "gene4"),
                translationDown = character(),
                translatedmRNAUp = character(),
                translatedmRNADown = character(),
                bufferingmRNAUp = character(),
                bufferingmRNADown = character(),
                mRNAAbundanceUp = character(),
                mRNAAbundanceDown = character(),
                totalmRNAUp = character(),
                totalmRNADown = character()
            ))
        },
        coloursSel = function(ads, genesIn, geneList, geneListcolours) "red",
        effectSel = function(ads, regulationGen, contrastSel, effectMeasure) c(gene1 = 1, gene2 = 2, gene3 = -1, gene4 = -2),
        getBg = function(ads, customBg, geneList) c("gene1", "gene2", "gene3", "gene4"),
        .package = "postNet"
    )

    ptn_load <- postNetStart(
        geneList = gene_list,
        geneListcolours = "blue",
        effectMeasure = effect_measure,
        source = "load",
        species = "human"
    )
    expect_equal(ptn_species(ptn_load), "human")
    expect_equal(ptn_version(ptn_load), "ver_40.202408")

    pos_file <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            id = c("tx1", "tx2", "tx3", "tx4"),
            UTR5_len = c(3, 3, 3, 3),
            CDS_stop = c(6, 6, 6, 6),
            Total_len = c(9, 9, 9, 9)
        ),
        file = pos_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    fasta_file <- tempfile(fileext = ".fa")
    writeLines(
        c(
            ">tx1", "AAATTTGGG",
            ">tx2", "CCCAAATTT",
            ">tx3", "GGGAAACCC",
            ">tx4", "TTTCCCGGG"
        ),
        fasta_file
    )

    testthat::local_mocked_bindings(
        read.fasta = function(file, seqtype = "AA") {
            list(
                tx1 = strsplit("AAATTTGGG", "")[[1]],
                tx2 = strsplit("CCCAAATTT", "")[[1]],
                tx3 = strsplit("GGGAAACCC", "")[[1]],
                tx4 = strsplit("TTTCCCGGG", "")[[1]]
            )
        },
        .package = "seqinr"
    )
    testthat::local_mocked_bindings(
        gffRead = function(gffFile) data.frame(dummy = 1),
        extGff = function(gff) data.frame(
            id = c("tx1", "tx2", "tx3", "tx4"),
            chr = "NC_1",
            strand = "+",
            start = 1,
            end = 9,
            transVer = 1,
            geneID = c("gene1", "gene2", "gene3", "gene4")
        ),
        gSel = function(annot, ads, customBg, geneList) annot,
        .package = "postNet"
    )

    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = "green",
            effectMeasure = effect_measure,
            source = "createFromFasta",
            posFile = pos_file,
            fastaFile = fasta_file,
            genomic_gff_file = tempfile(fileext = ".gff")
        ),
        "annotation are not compatible with gene IDs in the background"
    )

    src_gz1 <- tempfile(fileext = ".gbff.gz")
    src_gz2 <- tempfile(fileext = ".fa.gz")
    src_gz3 <- tempfile(fileext = ".gff.gz")
    file.create(src_gz1, src_gz2, src_gz3)

    testthat::local_mocked_bindings(
        gunzip = function(filename, remove = FALSE) sub("\\.gz$", "", filename),
        .package = "R.utils"
    )
    testthat::local_mocked_bindings(
        read.fasta = function(file, seqtype = "AA") {
            list(
                NM_001 = strsplit("AAATTTGGG", "")[[1]],
                NM_002 = strsplit("CCCAAATTT", "")[[1]],
                NM_003 = strsplit("GGGAAACCC", "")[[1]],
                NM_004 = strsplit("TTTCCCGGG", "")[[1]]
            )
        },
        .package = "seqinr"
    )
    testthat::local_mocked_bindings(
        system2 = function(command, args) {
            write.table(
                data.frame(
                    id = c("NM_001", "NM_002", "NM_003", "NM_004"),
                    UTR5_len = 3,
                    CDS_stop = 6,
                    Total_len = 9
                ),
                file = "customAnnot.txt",
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
            )
            0
        },
        file.rename = function(from, to) TRUE,
        .package = "base"
    )
    testthat::local_mocked_bindings(
        gffRead = function(gffFile) data.frame(dummy = 1),
        extGff = function(gff) data.frame(
            id = c("NM_001", "NM_002", "NM_003", "NM_004"),
            chr = "NC_1",
            strand = "+",
            start = 1,
            end = 9,
            transVer = 1,
            geneID = c("gene1", "gene2", "gene3", "gene4")
        ),
        .package = "postNet"
    )
    source_dir <- tempfile("sourcefiles")
    dir.create(source_dir)
    old_wd <- setwd(source_dir)
    on.exit(setwd(old_wd), add = TRUE)
    ptn_sourcefiles <- postNetStart(
        geneList = gene_list,
        geneListcolours = "orange",
        effectMeasure = effect_measure,
        source = "createFromSourceFiles",
        species = "human",
        rna_gbff_file = src_gz1,
        rna_fa_file = src_gz2,
        genomic_gff_file = src_gz3
    )
    expect_s4_class(ptn_sourcefiles, "postNetData")
})

test_that("postNetStart covers additional warnings and createFromFasta success paths", {
    gene_list <- list(groupA = c("gene1", "gene2"))
    effect_measure <- c(gene1 = 1, gene2 = -1)

    custom_file <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            id = c("tx1", "tx2"),
            geneID = c("gene1", "gene2"),
            UTR5_seq = c("AAA", "CCC"),
            CDS_seq = c("ATGAAA", "ATGCCC"),
            UTR3_seq = c("TTT", "AAA")
        ),
        file = custom_file,
        sep = "	",
        row.names = FALSE,
        quote = FALSE
    )

    expect_warning(
        ptn_warn <- postNetStart(
            geneList = gene_list,
            geneListcolours = "red",
            customBg = c("gene1", "gene2", "geneX"),
            effectMeasure = effect_measure,
            source = "custom",
            customFile = custom_file
        ),
        "genes in the background that are not present"
    )
    expect_s4_class(ptn_warn, "postNetData")

    pos_file <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            id = c("NM_001", "NM_002"),
            UTR5_len = c(3, 3),
            CDS_stop = c(6, 6),
            Total_len = c(9, 9)
        ),
        file = pos_file,
        sep = "	",
        row.names = FALSE,
        quote = FALSE
    )
    fasta_file <- tempfile(fileext = ".fa")
    writeLines(c(">NM_001", "AAATTTGGG", ">NM_002", "CCCAAATTT"), fasta_file)
    gff_file <- tempfile(fileext = ".gff")
    file.create(gff_file)

    testthat::local_mocked_bindings(
        read.fasta = function(file, seqtype = "AA") {
            list(
                NM_001 = strsplit("AAATTTGGG", "")[[1]],
                NM_002 = strsplit("CCCAAATTT", "")[[1]]
            )
        },
        .package = "seqinr"
    )
    testthat::local_mocked_bindings(
        gffRead = function(gffFile) data.frame(dummy = 1),
        extGff = function(gff) data.frame(
            id = c("NM_001", "NM_002"),
            chr = "NC_1",
            strand = "+",
            start = 1,
            end = 9,
            transVer = 1,
            geneID = c("gene1", "gene2")
        ),
        gSel = function(annot, ads, customBg, geneList) annot,
        .package = "postNet"
    )

    fasta_dir <- tempfile("createfasta")
    dir.create(fasta_dir)
    old <- setwd(fasta_dir)
    on.exit(setwd(old), add = TRUE)
    ptn_fasta <- postNetStart(
        geneList = gene_list,
        geneListcolours = "cyan",
        effectMeasure = effect_measure,
        source = "createFromFasta",
        posFile = pos_file,
        fastaFile = fasta_file,
        genomic_gff_file = gff_file,
        selection = "random"
    )
    expect_s4_class(ptn_fasta, "postNetData")
})

test_that("postNetStart covers remaining create, load, fasta, and adjustment branches", {
    gene_list <- list(groupA = c("gene1", "gene2"))
    effect_measure <- c(gene1 = 1, gene2 = -1)

    custom_file <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            id = c("tx1", "tx2"),
            geneID = c("gene1", "gene2"),
            UTR5_seq = c("AAA", "CCC"),
            CDS_seq = c("ATGAAA", "ATGCCC"),
            UTR3_seq = c("TTT", "AAA")
        ),
        file = custom_file,
        sep = "	",
        row.names = FALSE,
        quote = FALSE
    )

    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = "red",
            effectMeasure = effect_measure,
            source = "custom",
            customFile = custom_file,
            adjObj = list(UTR5 = c(tx1 = "GGG")),
            region_adj = "UTR5",
            excl = "yes"
        ),
        "input for 'excl'"
    )
    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = "red",
            effectMeasure = effect_measure,
            source = "custom",
            customFile = custom_file,
            adjObj = list(UTR5 = c(tx1 = "GGG")),
            region_adj = "UTR5",
            keepAll = "yes"
        ),
        "input for 'keepAll'"
    )

    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = "red",
            effectMeasure = effect_measure,
            source = "createFromSourceFiles",
            species = "rat",
            rna_gbff_file = tempfile(fileext = ".gbff.gz"),
            rna_fa_file = tempfile(fileext = ".fa.gz"),
            genomic_gff_file = tempfile(fileext = ".gff.gz")
        ),
        "Please specify a species"
    )

    testthat::local_mocked_bindings(
        getLink = function(url) NULL,
        .package = "postNet"
    )
    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = "red",
            effectMeasure = effect_measure,
            source = "create",
            species = "human"
        ),
        "current release directory"
    )

    testthat::local_mocked_bindings(
        getLink = function(url) if (grepl("current/$", url)) "GCF_000001405.40/" else NULL,
        .package = "postNet"
    )
    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = "red",
            effectMeasure = effect_measure,
            source = "create",
            species = "human"
        ),
        "version directory"
    )

    testthat::local_mocked_bindings(
        getLink = function(url) NULL,
        .package = "postNet"
    )
    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = "red",
            effectMeasure = effect_measure,
            source = "create",
            species = "mouse"
        ),
        "current release directory"
    )

    testthat::local_mocked_bindings(
        get_reference_data = function(file) data.frame(
            id = c("tx1", "tx2"),
            geneID = c("gene1", "gene2"),
            UTR5_seq = c("AAA", "CCC"),
            CDS_seq = c("ATGAAA", "ATGCCC"),
            UTR3_seq = c("TTT", "AAA")
        ),
        .package = "postNet"
    )
    ptn_load_mouse <- postNetStart(
        geneList = gene_list,
        geneListcolours = "red",
        effectMeasure = effect_measure,
        source = "load",
        species = "mouse"
    )
    expect_equal(ptn_version(ptn_load_mouse), "ver_27.202402")

    testthat::local_mocked_bindings(
        read.fasta = function(file, seqtype = "AA") {
            list(
                NM_001 = strsplit("AAATTTGGG", "")[[1]],
                NM_002 = strsplit("CCCAAATTT", "")[[1]]
            )
        },
        .package = "seqinr"
    )
    testthat::local_mocked_bindings(
        gffRead = function(gffFile) data.frame(dummy = 1),
        extGff = function(gff) data.frame(
            id = c("NM_001", "NM_002"),
            chr = "NC_1",
            strand = "+",
            start = 1,
            end = 9,
            transVer = 1,
            geneID = c("gene1", "gene2")
        ),
        gSel = function(annot, ads, customBg, geneList) annot,
        getLink = function(url) {
            if (grepl("annotation_releases/current/$", url)) {
                "GCF_000001635.27/"
            } else {
                "mock_genomic.gff.gz"
            }
        },
        download.file = function(url, destfile, ...) file.create(destfile),
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        download.file = function(url, destfile, ...) file.create(destfile),
        .package = "utils"
    )
    testthat::local_mocked_bindings(
        gunzip = function(filename, remove = TRUE) sub("\\.gz$", "", filename),
        .package = "R.utils"
    )
    testthat::local_mocked_bindings(
        file.remove = function(...) TRUE,
        .package = "base"
    )

    pos_file <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            id = c("NM_001", "NM_002"),
            UTR5_len = c(3, 3),
            CDS_stop = c(6, 6),
            Total_len = c(9, 9)
        ),
        file = pos_file,
        sep = "	",
        row.names = FALSE,
        quote = FALSE
    )
    fasta_file <- tempfile(fileext = ".fa")
    writeLines(c(">NM_001", "AAATTTGGG", ">NM_002", "CCCAAATTT"), fasta_file)

    fasta_dir <- tempfile("mousefasta")
    dir.create(fasta_dir)
    old <- setwd(fasta_dir)
    on.exit(setwd(old), add = TRUE)
    ptn_fasta_mouse <- postNetStart(
        geneList = gene_list,
        geneListcolours = "red",
        effectMeasure = effect_measure,
        source = "createFromFasta",
        posFile = pos_file,
        fastaFile = fasta_file,
        species = "mouse"
    )
    expect_s4_class(ptn_fasta_mouse, "postNetData")

    ptn_adjusted <- postNetStart(
        geneList = gene_list,
        geneListcolours = "red",
        effectMeasure = effect_measure,
        source = "custom",
        customFile = custom_file,
        adjObj = list(UTR5 = c(tx1 = "GGG")),
        region_adj = "UTR5",
        excl = FALSE,
        keepAll = TRUE
    )
    expect_equal(ptn_sequences(ptn_adjusted, region = "UTR5")[1], "GGG")
})

test_that("postNetStart validates inputs and can build from mocked create mode", {
    gene_list <- list(groupA = c("gene1", "gene2"))
    effect_measure <- c(gene1 = 1, gene2 = -1)

    expect_error(
        postNetStart(
            ads = mock_ads(),
            geneList = gene_list,
            source = "custom",
            customFile = tempfile(fileext = ".tsv")
        ),
        "either an anota2seq object or a gene list"
    )
    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = c("red", "blue"),
            effectMeasure = effect_measure,
            source = "custom",
            customFile = tempfile(fileext = ".tsv")
        ),
        "geneListcolours"
    )
    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = "red",
            customBg = 1:2,
            effectMeasure = effect_measure,
            source = "custom",
            customFile = tempfile(fileext = ".tsv")
        ),
        "customBg"
    )
    expect_error(
        postNetStart(
            geneList = gene_list,
            geneListcolours = "red",
            effectMeasure = effect_measure,
            source = "custom",
            customFile = tempfile(fileext = ".tsv"),
            adjObj = list(UTR5 = "AAA"),
            region_adj = "CDS"
        ),
        "region_adj"
    )

    testthat::local_mocked_bindings(
        getLink = function(url) {
            if (grepl("annotation_releases/current/$", url)) {
                "GCF_000001405.40/"
            } else {
                c("mock_rna.fna.gz", "mock_rna.gbff.gz", "mock_genomic.gff.gz")
            }
        },
        download.file = function(url, destfile, ...) file.create(destfile),
        gffRead = function(gffFile) data.frame(dummy = 1),
        extGff = function(gff) data.frame(
            id = c("NM_001", "NM_002"),
            chr = "NC_1",
            strand = "+",
            start = 1,
            end = 9,
            transVer = 1,
            geneID = c("gene1", "gene2")
        ),
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        download.file = function(url, destfile, ...) file.create(destfile),
        .package = "utils"
    )
    testthat::local_mocked_bindings(
        gunzip = function(filename, remove = TRUE) sub("\\.gz$", "", filename),
        .package = "R.utils"
    )
    testthat::local_mocked_bindings(
        read.fasta = function(file, seqtype = "AA") {
            list(
                NM_001 = strsplit("AAATTTGGG", "")[[1]],
                NM_002 = strsplit("CCCAAATTT", "")[[1]]
            )
        },
        .package = "seqinr"
    )
    testthat::local_mocked_bindings(
        system2 = function(command, args) {
            write.table(
                data.frame(
                    id = c("NM_001", "NM_002"),
                    UTR5_len = 3,
                    CDS_stop = 6,
                    Total_len = 9
                ),
                file = "customAnnot.txt",
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
            )
            0
        },
        .package = "base"
    )

    create_dir <- tempfile("create")
    dir.create(create_dir)
    old_wd <- setwd(create_dir)
    on.exit(setwd(old_wd), add = TRUE)
    ptn_create <- postNetStart(
        geneList = gene_list,
        geneListcolours = "purple",
        effectMeasure = effect_measure,
        source = "create",
        species = "human"
    )
    expect_s4_class(ptn_create, "postNetData")
    expect_equal(ptn_version(ptn_create), "GCF_000001405.40")
})

test_that("postNetStart handles ads-based inputs and mouse create mode", {
    ads <- mock_ads_plot()
    effect_df <- data.frame(apvEff = c(1, -1, 0.5, -0.5), row.names = c("gene1", "gene2", "gene3", "gene4"))
    if ("totalmRNA" %in% slotNames(ads)) {
        ads@totalmRNA@apvStatsRvm <- list(effect_df)
        ads@translatedmRNA@apvStatsRvm <- list(effect_df)
        ads@translation@apvStatsRvm <- list(effect_df)
        ads@buffering@apvStatsRvm <- list(effect_df)
    }

    expect_error(
        postNetStart(
            ads = ads,
            regulation = "badMode",
            contrast = 1,
            source = "load",
            species = "human"
        ),
        "regulation"
    )
    expect_error(
        postNetStart(
            ads = ads,
            regulation = "translationUp",
            contrast = 3,
            source = "load",
            species = "human"
        ),
        "The input for 'contrast' must be a numeric vector"
    )

    testthat::local_mocked_bindings(
        getLink = function(url) {
            if (grepl("annotation_releases/current/$", url)) {
                "GCF_000001635.27/"
            } else {
                c("mock_rna.fna.gz", "mock_rna.gbff.gz", "mock_genomic.gff.gz")
            }
        },
        download.file = function(url, destfile, ...) file.create(destfile),
        gffRead = function(gffFile) data.frame(dummy = 1),
        extGff = function(gff) data.frame(
            id = c("NM_001", "NM_002"),
            chr = "NC_1",
            strand = "+",
            start = 1,
            end = 9,
            transVer = 1,
            geneID = c("gene1", "gene2")
        ),
        .package = "postNet"
    )
    testthat::local_mocked_bindings(
        gunzip = function(filename, remove = TRUE) sub("\\.gz$", "", filename),
        .package = "R.utils"
    )
    testthat::local_mocked_bindings(
        read.fasta = function(file, seqtype = "AA") {
            list(
                NM_001 = strsplit("AAATTTGGG", "")[[1]],
                NM_002 = strsplit("CCCAAATTT", "")[[1]]
            )
        },
        .package = "seqinr"
    )
    testthat::local_mocked_bindings(
        system2 = function(command, args) {
            write.table(
                data.frame(
                    id = c("NM_001", "NM_002"),
                    UTR5_len = 3,
                    CDS_stop = 6,
                    Total_len = 9
                ),
                file = "customAnnot.txt",
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
            )
            0
        },
        .package = "base"
    )
    create_dir <- tempfile("createmouse")
    dir.create(create_dir)
    old_wd <- setwd(create_dir)
    on.exit(setwd(old_wd), add = TRUE)
    ptn_mouse <- postNetStart(
        geneList = list(groupA = c("gene1", "gene2")),
        geneListcolours = "brown",
        effectMeasure = c(gene1 = 1, gene2 = -1),
        source = "create",
        species = "mouse"
    )
    expect_s4_class(ptn_mouse, "postNetData")
    expect_equal(ptn_version(ptn_mouse), "GCF_000001635.27")
})
