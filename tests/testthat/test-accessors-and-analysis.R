test_that("core postNet accessors return stored slots", {
    ptn <- mock_ptn(
        with_codon = TRUE,
        with_miRNA = TRUE,
        with_gsea = TRUE,
        with_gage = TRUE,
        with_motifs = TRUE
    )

    expect_equal(ptn_sequences(ptn, "UTR3"), c("TTT", "AAA"))
    expect_equal(ptn_id(ptn, "CDS"), c("tx1", "tx2"))
    expect_equal(ptn_geneID(ptn, "UTR5"), c("gene1", "gene2"))
    expect_equal(ptn_geneList(ptn), list(comparisonA = c("gene1", "gene2")))
    expect_equal(ptn_background(ptn), c("gene1", "gene2", "gene3"))
    expect_equal(ptn_effect(ptn), c(gene1 = 1.2, gene2 = -0.4))
    expect_equal(ptn_colours(ptn), "red")
    expect_equal(ptn_species(ptn), "human")
    expect_equal(ptn_version(ptn), "v1")
    expect_equal(ptn_selection(ptn), "random")
    expect_equal(
        ptn_motifSelection(ptn, "UTR3"),
        data.frame(score = 2, row.names = "motif2")
    )
    expect_equal(
        ptn_motifGeneList(ptn, "UTR3", "comparisonA"),
        data.frame(score = 2, row.names = "motif2")
    )
    expect_equal(ptn_codonSelection(ptn, "comparisonA"), c("AAA", "AAG"))
    expect_equal(
        ptn_features(ptn),
        data.frame(featureA = c(1, 2), row.names = c("gene1", "gene2"))
    )
})

test_that("core accessors validate their inputs", {
    ptn <- mock_ptn()
    ptn_with_motifs <- mock_ptn(with_motifs = TRUE)

    expect_error(ptn_sequences(ptn, c("UTR3", "CDS")), "can only be one of")
    expect_error(
        ptn_motifGeneList(ptn_with_motifs, "UTR3", "missingComparison"),
        "is not stored in the postNetData object"
    )
    expect_error(ptn_geneID(list(), "UTR3"), "unable to find an inherited method")
})

test_that("signCalc builds named binary signatures", {
    ptn <- mock_ptn()
    signatures <- list(sigA = "gene1", sigB = c("gene2", "gene3"))

    out <- signCalc(ptn, signatures)

    expect_named(out, c("sigA", "sigB"))
    expect_equal(out$sigA, c(gene1 = 1, gene2 = 0))
    expect_equal(out$sigB, c(gene1 = 0, gene2 = 1))
})

test_that("analysis accessors filter and reshape results", {
    ptn <- mock_ptn(
        with_codon = TRUE,
        with_miRNA = TRUE,
        with_gsea = TRUE,
        with_gage = TRUE
    )

    codon_df <- ptn_codonAnalysis(ptn)
    expect_s3_class(codon_df, "data.frame")
    expect_equal(
        colnames(codon_df),
        c(
            "geneID",
            "codon",
            "AA",
            "count",
            "frequency",
            "AACountPerGene",
            "relative_frequency"
        )
    )
    expect_equal(nrow(codon_df), 3)

    mirna_res <- ptn_miRNA_analysis(ptn, direction = "greater", threshold = 0.05)
    expect_equal(mirna_res$id, "mir1")
    expect_equal(
        colnames(mirna_res),
        c("id", "score", "context", "qvalue", "support", "pvalue")
    )

    expect_equal(
        ptn_miRNA_to_gene(ptn, c("mir1", "missing")),
        list(mir1 = c("gene1", "gene2"))
    )

    gsea_res <- ptn_GSEA(ptn, threshold = 0.05)
    expect_equal(gsea_res$pathway, "setA")

    gage_res <- ptn_GAGE(ptn, category = "BP", direction = "greater", threshold = 0.05)
    expect_equal(rownames(gage_res), "termA")
})

test_that("analysis accessors report missing analyses and invalid values", {
    ptn <- mock_ptn()

    expect_error(
        ptn_miRNA_analysis(ptn, direction = "greater", threshold = 0.05),
        "Please run the miRNAanalysis\\(\\) function first"
    )
    expect_error(ptn_GSEA(ptn), "Please run the gseaAnalysis\\(\\) function first")
    expect_error(
        ptn_GAGE(ptn, category = "BP", direction = "greater", threshold = 0.05),
        "Please run the gageAnalysis\\(\\) function first"
    )
    expect_error(
        ptn_miRNA_analysis(
            mock_ptn(with_miRNA = TRUE),
            direction = "up",
            threshold = 0.05
        ),
        "must be either 'greater' or 'less'"
    )
})

test_that("postNet analysis accessors cover GO, empty outputs and model helpers", {
    ptn <- mock_ptn(
        with_miRNA = TRUE,
        with_gsea = TRUE,
        with_gage = TRUE,
        with_feature_models = TRUE
    )
    ptn@analysis@GO <- methods::new(
        "postNetGO",
        BP = list(
            comparisonA = methods::new(
                "mockEnrichResult",
                result = data.frame(
                    ID = c("GO:1", "GO:2"),
                    Description = c("term1", "term2"),
                    p.adjust = c(0.01, 0.2),
                    geneID = c("gene1:gene2", "gene2"),
                    row.names = c("row1", "row2")
                )
            )
        ),
        CC = NULL,
        MF = NULL,
        KEGG = NULL
    )

    go_res <- ptn_GO(ptn, category = "BP", geneList = "comparisonA", threshold = 0.05)
    expect_s3_class(go_res, "data.frame")
    expect_equal(go_res$ID, "GO:1")
    expect_message(
        ptn_GO(ptn, category = "BP", geneList = "comparisonA", threshold = 1e-10),
        "There are no enriched GO categories to output"
    )

    expect_message(
        ptn_miRNA_analysis(ptn, direction = "less", threshold = 1e-10),
        "There are no enriched miRNAs to output"
    )
    expect_message(
        ptn_GSEA(ptn, threshold = 1e-10),
        "There are no enriched gene sets to output from GSEA"
    )
    expect_message(
        ptn_GAGE(ptn, category = "BP", direction = "less", threshold = 0.05),
        "There are no enriched terms to output from GAGE"
    )

    expect_output(ptn_check_models(ptn, "lm"), "NULL")
    expect_error(ptn_check_models(mock_ptn(), "lm"), "Please run lm analysis first")

    expect_error(ptn_GO(ptn, category = c("BP", "CC"), geneList = "comparisonA", threshold = 0.05), "Please provide only one GO category")
    expect_error(ptn_GO(ptn, category = "BP", geneList = "missing", threshold = 0.05), "included in the postNetData object")
    expect_error(ptn_GO(ptn, category = "BP", geneList = "comparisonA", threshold = "low"), "single numeric value")

    expect_error(ptn_GSEA(ptn, threshold = "low"), "single numeric value")
    expect_error(ptn_GAGE(ptn, category = c("BP", "CC"), direction = "greater", threshold = 0.05), "Please provide only one category")
    expect_error(ptn_GAGE(ptn, category = "BP", direction = "greater", threshold = "low"), "single numeric value")
    expect_error(ptn_model(ptn, analysis_type = "rf", model = "finalModel", comparison = "one"), "single numeric value")
    expect_error(ptn_selectedFeatures(ptn, analysis_type = "rf", comparison = "one"), "single numeric value")
    expect_error(ptn_networkGraph(ptn, comparison = "one"), "single numeric value")
})
