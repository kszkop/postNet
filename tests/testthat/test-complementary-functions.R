test_that("GFF and sequence helpers process compact examples", {
    attrs <- c(
        "gene=GENE1;transcript_id=NM_000001.2",
        "gene=GENE2;transcript_id=XM_000002.1"
    )
    expect_equal(postNet:::getAttributeField(attrs, "gene"), c("GENE1", "GENE2"))
    expect_equal(
        postNet:::getAttributeField(attrs, "missing"),
        c(NA_character_, NA_character_)
    )

    gff <- data.frame(
        seqname = c("NC_1", "NC_2", "NC_3"),
        source = "src",
        feature = c("mRNA", "transcript", "exon"),
        start = c(1L, 5L, 10L),
        end = c(9L, 12L, 13L),
        score = ".",
        strand = c("+", "-", "+"),
        frame = ".",
        attributes = c(
            "gene=GENE1;transcript_id=NM_000001.2",
            "gene=GENE2;transcript_id=XM_000002.1",
            "gene=GENE3;transcript_id=NM_000003.1"
        )
    )

    bed <- postNet:::extGff(gff)
    expect_equal(colnames(bed), c("id", "chr", "strand", "start", "end", "transVer", "geneID"))
    expect_equal(bed$id, "NM_000001")
    expect_equal(bed$geneID, "GENE1")

    expect_equal(postNet:::extract_seq(0, "ATGAAA"), "ATG")
    expect_equal(
        postNet:::generate_pairs(c(0, 1, 2)),
        stats::setNames(list(c(0, 1), c(0, 2), c(1, 2)), c("1", "2", "3"))
    )
})

test_that("annotation selection helpers subset, measure, and adjust transcripts", {
    annot <- data.frame(
        id = c("tx1", "tx2", "tx3"),
        geneID = c("gene1", "gene1", "gene2"),
        UTR5_seq = c("AA", "AAAA", "A"),
        CDS_seq = c("ATG", "ATGAAA", "AT"),
        UTR3_seq = c("TT", "TTT", "T"),
        stringsAsFactors = FALSE
    )

    geneList <- list(groupA = c("gene1", "gene2"))

    expect_equal(
        postNet:::gSel(annot, ads = NULL, customBg = c("gene2"), geneList = geneList)$geneID,
        c("gene1", "gene1", "gene2")
    )
    expect_equal(
        postNet:::gSel(annot, ads = NULL, customBg = NULL, geneList = list(groupA = "gene2"))$geneID,
        "gene2"
    )

    reg <- postNet:::regSel(annot, "CDS")
    expect_equal(colnames(reg), c("id", "geneID", "seqTmp", "lenTmp"))
    expect_equal(reg$lenTmp, c(3, 6, 2))

    expect_equal(postNet:::isoSel(reg, "shortest")$id, c("tx1", "tx3"))
    expect_equal(postNet:::isoSel(reg, "longest")$id, c("tx2", "tx3"))
    expect_equal(
        postNet:::isoSel(reg, "random", setSeed = 99)$id,
        postNet:::isoSel(reg, "random", setSeed = 99)$id
    )

    adjusted <- postNet:::adjustSeq(
        annot,
        adjObj = list(UTR3 = c(tx2 = "GGG")),
        region_adj = "UTR3",
        excl = FALSE,
        keepAll = FALSE
    )
    expect_equal(adjusted$UTR3_seq[adjusted$id == "tx2"], "GGG")
    expect_false("tx1" %in% adjusted$id)

    adjusted_excl <- postNet:::adjustSeq(
        annot,
        adjObj = list(UTR5 = c(tx3 = "CCCC")),
        region_adj = "UTR5",
        excl = TRUE,
        keepAll = TRUE
    )
    expect_equal(adjusted_excl$id, "tx3")
    expect_equal(adjusted_excl$UTR5_seq, "CCCC")
})

test_that("background, feature prep, and identifier helpers are deterministic", {
    ptn <- mock_ptn(with_background = TRUE)

    expect_equal(
        postNet:::getBg(
            ads = NULL,
            customBg = c("gene1", "gene2"),
            geneList = list(groupA = "gene1")
        ),
        c("gene1", "gene2")
    )
    expect_error(
        postNet:::getBg(
            ads = NULL,
            customBg = c("gene2"),
            geneList = list(groupA = "gene1")
        ),
        "not present in 'customBg'"
    )
    expect_equal(
        postNet:::getBg(
            ads = NULL,
            customBg = NULL,
            geneList = list(groupA = c("gene1", "gene2"))
        ),
        c("gene1", "gene2")
    )

    expect_equal(postNet:::check_id_type("1234"), "entrezID")
    expect_equal(postNet:::check_id_type("TP53"), "geneID")
    expect_equal(postNet:::check_id_type("ENSG-1"), "unknown")

    prepped <- suppressMessages(postNet:::prepFeatures(
        ptn,
        list(
            featureA = c(gene1 = 1, gene2 = 2),
            featureB = c(gene1 = NA, gene2 = 3)
        )
    ))
    expect_equal(rownames(ptn_features(prepped)), "gene2")
    expect_equal(colnames(ptn_features(prepped)), c("featureA", "featureB"))
})


test_that("ads helpers and parser utilities in complementaryFunctions work on mocked inputs", {
    if (!methods::isClass("MockSelectedRvm")) {
        methods::setClass("MockSelectedRvm", slots = c(selectedRvmData = "list"))
    }
    if (!methods::isClass("MockAbundanceRvm")) {
        methods::setClass("MockAbundanceRvm", slots = c(totalmRNA = "list"))
    }
    if (!methods::isClass("MockFullAds")) {
        methods::setClass(
            "MockFullAds",
            slots = c(
                contrasts = "matrix",
                dataP = "matrix",
                totalmRNA = "ANY",
                translatedmRNA = "ANY",
                buffering = "ANY",
                translation = "ANY",
                selectedTranslation = "MockSelectedRvm",
                selectedTranslatedmRNA = "MockSelectedRvm",
                selectedBuffering = "MockSelectedRvm",
                mRNAAbundance = "MockAbundanceRvm",
                selectedTotalmRNA = "MockSelectedRvm"
            )
        )
    }

    reg_df <- data.frame(
        apvEff = c(1, 2, 3, 4, -1, -2, -3, -4),
        singleRegMode = rep("translation", 8),
        row.names = paste0("gene", 1:8)
    )
    total_df <- data.frame(
        apvEff = c(1, 2, 3, 4, -1, -2, -3, -4),
        singleRegMode = rep("abundance", 8),
        row.names = paste0("gene", 9:16)
    )
    ads_full <- methods::new(
        "MockFullAds",
        contrasts = matrix(1, nrow = 1, ncol = 2),
        dataP = matrix(1:16, ncol = 1, dimnames = list(paste0("gene", 1:16), "x")),
        totalmRNA = NULL,
        translatedmRNA = NULL,
        buffering = NULL,
        translation = NULL,
        selectedTranslation = methods::new("MockSelectedRvm", selectedRvmData = list(reg_df, reg_df)),
        selectedTranslatedmRNA = methods::new("MockSelectedRvm", selectedRvmData = list(reg_df, reg_df)),
        selectedBuffering = methods::new("MockSelectedRvm", selectedRvmData = list(reg_df, reg_df)),
        mRNAAbundance = methods::new("MockAbundanceRvm", totalmRNA = list(total_df, total_df)),
        selectedTotalmRNA = methods::new("MockSelectedRvm", selectedRvmData = list(total_df, total_df))
    )

    regs <- postNet:::anota2seqGetDirectedRegulations(ads_full)
    expect_equal(length(regs), 2)
    expect_true("gene1" %in% regs[[1]]$translationUp)
    expect_true("gene5" %in% regs[[1]]$translationDown)
    expect_true("gene9" %in% regs[[1]]$totalmRNAUp)
    expect_true("gene13" %in% regs[[1]]$totalmRNADown)

    res_filtered <- postNet:::resSel(
        ads = ads_full,
        regulation = c("translationUp", "translationDown"),
        contrast = c(1, 2),
        geneList = NULL
    )
    expect_named(res_filtered, c("translationUp_c1", "translationDown_c2"))

    res_all <- postNet:::resSel(ads = ads_full, regulation = NULL, contrast = NULL, geneList = NULL)
    expect_true(any(grepl("translationUp_c1", names(res_all))))

    cols_ads <- postNet:::coloursSel(ads_full, genesIn = res_filtered, geneList = NULL, geneListcolours = NULL)
    expect_length(cols_ads, 2)

    testthat::local_mocked_bindings(
        anota2seqGetOutput = function(ads, output, selContrast, getRVM) {
            data.frame(
                identifier = c("gene1", "gene2"),
                translation.apvEff = c(1.1, -1.1),
                totalmRNA.apvEff = c(0.5, -0.5)
            )
        },
        .package = "anota2seq"
    )
    eff_translation <- postNet:::effectSel(ads_full, regulationGen = "translation", contrastSel = 1, effectMeasure = NULL)
    expect_equal(unname(eff_translation), c(1.1, -1.1))
    eff_abundance <- postNet:::effectSel(ads_full, regulationGen = "mRNAAbundance", contrastSel = 1, effectMeasure = NULL)
    expect_equal(unname(eff_abundance), c(0.5, -0.5))

    gff_file <- tempfile(fileext = ".gff")
    writeLines(
        c(
            "NC_1	src	mRNA	1	9	.	+	.	gene=GENE1;transcript_id=NM_000001.2",
            "NC_1	src	transcript	5	12	.	+	.	gene=GENE2;transcript_id=NM_000002.1"
        ),
        gff_file
    )
    expect_message(gff_parsed <- postNet:::gffRead(gff_file), "Reading")
    expect_equal(nrow(gff_parsed), 2)
})
