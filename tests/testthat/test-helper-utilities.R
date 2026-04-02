test_that("utility helpers in complementaryFunctions behave deterministically", {
    expect_equal(postNet:::antilog(log2(8), 2), 8)
    expect_true(postNet:::is_by_3(c("ATG", "AAATTT")))
    expect_false(postNet:::is_by_3(c("AT", "AAA")))

    expect_equal(postNet:::adjust_ylim(-2, 3), c(-2, 3))
    expect_equal(postNet:::adjust_ylim(-5, -1), c(-5, 0))
    expect_equal(postNet:::adjust_ylim(1, 4), c(0, 4))

    expect_equal(postNet:::motifLenCalc("A[CG]T"), 3)
    expect_equal(postNet:::replaceProtAmbig("AX"), "A[ACDEFGHIKLMNPQRSTVWY]")
    expect_equal(unname(postNet:::remove_last3(c("ATGAAA", "CCCTTT"))), c("ATG", "CCC"))

    cc1 <- postNet:::codonCount("ATGAAA", "gene1", codonN = 1)
    expect_true(all(c("geneID", "codon", "AA", "count", "frequency", "AACountPerGene", "relative_frequency") %in% colnames(cc1)))
    expect_true(any(cc1$codon == "ATG"))

    cc2 <- postNet:::codonCount("ATGAAA", "gene1", codonN = 2)
    expect_true(all(c("geneID", "codon", "count", "AA", "frequency") %in% colnames(cc2)))

    codon_df <- data.frame(
        codon = c("AAA", "AAG", "GGA", "GGG"),
        AA = c("K", "K", "G", "G"),
        A = c(10, 5, 4, 2),
        B = c(2, 8, 1, 5)
    )
    aa_df <- data.frame(
        AA = c("K", "G"),
        A = c(15, 6),
        B = c(10, 6)
    )
    expect_true(all(names(postNet:::statOnDf(codon_df, c("A", "B"), "codon")) %in% c("AAA", "AAG", "GGA", "GGG")))
    expect_true(all(names(postNet:::statOnDf(aa_df, c("A", "B"), "AA")) %in% c("K", "G")))

    expect_equal(postNet:::roundNice(3.2, direction = "up"), 4)
    expect_equal(postNet:::roundNice(3.2, direction = "down"), 2)
    expect_error(postNet:::roundNice(c(1, 2), direction = "up"), "'x' must be of length 1")

    expect_equal(postNet:::combSeq(list(c("AAA", "CCC"), c("TTT", "GGG"))), list("AAATTT", "CCCGGG"))

    uorf_n <- postNet:::calc_uORF("CCCCCCC", "TAA", "[ATGC][ATGC][ATGC]ATG[ATGC]", "number")
    expect_equal(uorf_n, 0)
    expect_true(is.na(postNet:::calc_uORF("CCCCCCC", "TAA", "[ATGC][ATGC][ATGC]ATG[ATGC]", "position")))

    expect_equal(postNet:::dup_tolerance(c(1, 1 + 1e-9, 3)), c(TRUE, TRUE, FALSE))
    dups <- postNet:::getDup(stats::setNames(c(1, 1 + 1e-9, 3), c("a", "b", "c")))
    expect_length(dups, 1)
    expect_message(postNet:::printDup(dups), "too highly corelated")

    expect_equal(postNet:::rescale(5, 0, 10, 0, 1), 0.5)
    expect_match(postNet:::wrapNames("a very long name", 4), "\n")

    g <- igraph::make_ring(3)
    lay <- postNet:::layoutCalc(g, n = 4)
    expect_equal(dim(lay), c(3, 2))

    annot_seq <- data.frame(
        seq = c("AAATTTGGG", "CCCAAATTT"),
        UTR5_len = c(3, 3),
        CDS_stop = c(6, 6),
        Total_len = c(9, 9)
    )
    reg_seqs <- postNet:::extractRegSeq(annot_seq)
    expect_equal(reg_seqs$UTR5_seq, c("AAA", "CCC"))
    expect_equal(reg_seqs$CDS_seq, c("TTT", "AAA"))
    expect_equal(reg_seqs$UTR3_seq, c("GGG", "TTT"))

    expect_match(postNet:::generateOut("id1", list(id1 = c(a = 1, b = 2))), "id1")

    layout <- matrix(c(0, 0, 10, 20), ncol = 2, byrow = TRUE)
    norm_layout <- postNet:::normalizeLayout(layout)
    expect_equal(range(norm_layout[, 1]), c(-1, 1))
    expect_equal(range(norm_layout[, 2]), c(-1, 1))

    expect_equal(unname(postNet:::colourAssign(c("A", "B"), c(A = "red", B = "blue"))), c("red", "blue"))
    expect_true(postNet:::is_binary(c(0, 1, 0, 1)))
    expect_false(postNet:::is_binary(c(0, 2)))
    expect_true(postNet:::is_categorical(c(1, 1, 2, 2)))
    expect_false(postNet:::is_categorical(runif(20)))
})

test_that("plot_fmap and cache/web helpers can be exercised with mocks", {
    fmap <- data.frame(UMAP1 = c(0, 1), UMAP2 = c(1, 0))

    pf_num <- postNet:::plot_fmap(fmap, c(-1, 2), name = "Effect")
    expect_true(all(c("mainPlot", "legend") %in% names(pf_num)))

    pf_bin <- postNet:::plot_fmap(fmap, c(0, 1), name = "Binary")
    expect_true(all(c("mainPlot", "legend") %in% names(pf_bin)))

    local({
        testthat::local_mocked_bindings(
            request = function(url) structure(list(url = url), class = "mockreq"),
            req_perform = function(req) "response",
            resp_check_status = function(resp) resp,
            resp_body_string = function(resp) "<html><body><a href='a.txt'>a</a><a href='b.txt'>b</a></body></html>",
            .package = "httr2"
        )
        expect_equal(postNet:::getLink("https://example.org"), c("a.txt", "b.txt"))
    })

    local({
        testthat::local_mocked_bindings(
            request = function(url) stop("boom"),
            .package = "httr2"
        )
        expect_warning(expect_null(postNet:::getLink("https://example.org")), "Failed to fetch")
    })

    tmp_gz <- tempfile(fileext = ".txt.gz")
    con <- gzfile(tmp_gz, "wt")
    writeLines("id\tgeneID\tCDS_seq\nx\tg\tATG", con)
    close(con)

    testthat::local_mocked_bindings(
        BiocFileCache = function(cache_dir, ask = FALSE) "bfc",
        bfcquery = function(bfc, query, field = "rname", exact = TRUE) data.frame(rid = "RID1"),
        bfcrpath = function(bfc, rids = NULL, rid = NULL) tmp_gz,
        .package = "BiocFileCache"
    )
    expect_equal(postNet:::get_reference_data("file.txt.gz")$id, "x")

    removed <- character()
    testthat::local_mocked_bindings(
        BiocFileCache = function(cache_dir, ask = FALSE) "bfc",
        bfcinfo = function(bfc) data.frame(rid = c("a", "b")),
        bfcremove = function(bfc, rid) {
            removed <<- rid
            invisible(TRUE)
        },
        .package = "BiocFileCache"
    )
    expect_true(postNet:::clear_postNet_cache(tempdir()))
    expect_equal(removed, c("a", "b"))
})

test_that("plotting and modeling helpers can run on compact synthetic inputs", {
    out_dir <- tempfile("helpers")
    dir.create(out_dir)
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)

    pdf("plotpostnet_box.pdf")
    postNet:::plotPostNet(
        resOut = list(background = c(1, 2, 3, 4), comparisonA = c(2, 3, 4, 5)),
        colOut = c("grey45", "red"),
        comparisons = list(c(0, 1)),
        ylabel = "test",
        plotType = "boxplot"
    )
    dev.off()
    expect_true(file.exists(file.path(out_dir, "plotpostnet_box.pdf")))

    pdf("plotpostnet_ecdf.pdf")
    postNet:::plotPostNet(
        resOut = list(background = c(1, 2, 3, 4), comparisonA = c(2, 3, 4, 5)),
        colOut = c("grey45", "red"),
        comparisons = list(c(0, 1)),
        ylabel = "test",
        plotType = "ecdf"
    )
    dev.off()
    expect_true(file.exists(file.path(out_dir, "plotpostnet_ecdf.pdf")))

    postNet:::plotScatterInd(
        set1 = data.frame(x = c(1, 2, 3), y = c(2, 3, 5)),
        set2 = data.frame(x = c(2, 4, 5), y = c(1, 3, 4)),
        orgName = "feat",
        coloursIn = c("red", "blue"),
        nameOut = "scatter"
    )
    expect_true(file.exists(file.path(out_dir, "scatter_feat_individually.pdf")))

    dataIn <- data.frame(
        a1 = c(1, 2, 3, 4, 5, 6),
        a2 = c(1, 1, 2, 3, 5, 8),
        effM = c(1, 2, 1.5, 3, 2.5, 4)
    )
    namesDf <- data.frame(
        originalNames = c("feat1", "feat2"),
        newNames = c("a1", "a2")
    )

    testthat::local_mocked_bindings(
        plot.igraph = function(...) NULL,
        .package = "igraph"
    )
    lm_out <- suppressMessages(postNet:::runLM(
        dataIn = dataIn,
        namesDf = namesDf,
        allFeat = TRUE,
        useCorel = FALSE,
        covarFilt = 20,
        nameOut = "runlm",
        NetModelSel = "adjusted",
        coloursIn = c("red", "blue"),
        lmfeatGroup = NULL,
        lmfeatGroupColour = NULL,
        fdrUni = 1,
        stepP = 1
    ))
    expect_s4_class(lm_out, "postNetFeatureIntegration_lm")
    expect_true(file.exists(file.path(out_dir, "runlm_FinalModel.pdf")))
    expect_true(file.exists(file.path(out_dir, "runlm_network.pdf")))
})

test_that("signature and identifier helpers reach real package data", {
    human_sigs <- postNet:::get_signatures("human")
    mouse_sigs <- postNet:::get_signatures("mouse")
    expect_true(is.list(human_sigs))
    expect_true(is.list(mouse_sigs))
    expect_error(postNet:::get_signatures("rat"), "Currently, 'human' or 'mouse' are available")

    expect_true(postNet:::check_featSel(c("a", "b"), data.frame(a = 1, b = 2)))
    expect_false(postNet:::check_featSel("a", data.frame(a = 1, b = 2)))
    expect_true(postNet:::check_featCol("a", data.frame(a = 1, b = 2)))
    expect_true(postNet:::check_predFeat(data.frame(a = c(1, 2), row.names = c("g1", "g2"))))
    expect_false(postNet:::check_predFeat(data.frame(a = c(1, NA), row.names = c("g1", "g2"))))
})
