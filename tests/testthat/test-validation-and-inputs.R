test_that("validation helpers accept valid inputs and reject invalid ones", {
    expect_null(postNet:::check_region("UTR3"))
    expect_error(postNet:::check_region("bad"), "must contain valid values")
    expect_error(postNet:::check_region(NULL), "must be a non-empty character")

    expect_null(postNet:::check_adjObj(list(UTR5 = "AAA", UTR3 = "TTT")))
    expect_error(postNet:::check_adjObj("bad"), "not a list")
    expect_error(postNet:::check_adjObj(list(CDS = "ATG")), "only be 'UTR3' or 'UTR5'")

    expect_null(postNet:::check_selection("Longest"))
    expect_error(postNet:::check_selection("middle"), "must be one of")

    expect_null(postNet:::check_plotType("Violin"))
    expect_error(postNet:::check_plotType(NULL), "must be one of")

    expect_null(postNet:::checkAnnot(data.frame(
        id = "tx1",
        geneID = "gene1",
        UTR5_seq = "AAA",
        CDS_seq = "ATG",
        UTR3_seq = "TTT"
    )))
    expect_error(postNet:::checkAnnot(data.frame(id = "tx1")), "missing in 'annot'")

    expect_null(postNet:::checkAnnotCod(data.frame(
        id = "tx1",
        geneID = "gene1",
        CDS_seq = "ATG"
    )))
    expect_error(
        postNet:::checkAnnotCod(data.frame(id = "tx1")),
        "missing in 'customFileCod'"
    )

    expect_true(postNet:::check_comparisons(list(c(0, 1), c(1, 2))))
    expect_false(postNet:::check_comparisons("bad"))

    expect_true(postNet:::is_valid_named_list(list(a = c(1, 2), b = 3)))
    expect_false(postNet:::is_valid_named_list(list(c(1, 2))))

    expect_true(postNet:::is_numeric_vector(c(1, 2)))
    expect_false(postNet:::is_numeric_vector(list(1, 2)))

    expect_true(postNet:::check_logical(TRUE))
    expect_error(postNet:::check_logical(c(TRUE, FALSE)), "length = 2")
    expect_true(postNet:::check_number(1.5))
    expect_false(postNet:::check_number(NA_real_))

    expect_true(postNet:::is_named_list_of_named_numeric_vectors(
        list(a = c(x = 1, y = 2))
    ))
    expect_false(postNet:::is_named_list_of_named_numeric_vectors(list(a = c(1, 2))))

    expect_null(postNet:::check_source("load"))
    expect_error(postNet:::check_source("remote"), "Invalid source")
    expect_null(postNet:::checkSourceFE("custom"))
    expect_error(postNet:::checkSourceFE("loadFromWeb"), "Invalid 'sourceFE'")

    expect_true(postNet:::is_valid_species("Human"))
    expect_false(postNet:::is_valid_species("rat"))
    expect_null(postNet:::check_DNAsequence(c("ACGT", "ACGT123")))
    expect_error(postNet:::check_DNAsequence(1), "must be a named list of")
    expect_error(
        postNet:::check_DNAsequence(c("ACGU")),
        "do not appear to all be DNA sequences"
    )
    expect_true(postNet:::is_valid_seq_type("RNA"))
    expect_false(postNet:::is_valid_seq_type("lipid"))
    expect_true(postNet:::isStartCodon("ATG"))
    expect_false(postNet:::isStartCodon("AT"))
})

test_that("codon and analysis validation helpers behave as expected", {
    expect_true(postNet:::is_valid_analysis("codon"))
    expect_false(postNet:::is_valid_analysis("protein"))
    expect_true(postNet:::isUnit("count"))
    expect_false(postNet:::isUnit("counts"))

    valid_codon_df <- data.frame(
        geneID = "gene1",
        codon = "AAA",
        AA = "Lys",
        count = 1,
        frequency = 0.1,
        AACountPerGene = 1,
        relative_frequency = 0.1
    )

    expect_null(postNet:::check_codonIn(valid_codon_df))
    expect_error(
        postNet:::check_codonIn(valid_codon_df[, 1:3]),
        "does not contain all required elements"
    )
    expect_true(postNet:::check_codons(list(lys = c("AAA", "AAG"))))
    expect_false(postNet:::check_codons(list(lys = "XYZ")))
    expect_true(postNet:::check_AA(list(lys = c("K", "Lys"))))
    expect_false(postNet:::check_AA(list(lys = "Xxx")))
    expect_true(postNet:::checkSlopes(-1, 2))
    expect_error(postNet:::checkSlopes(NA_real_, 2), "minSlope")
    expect_null(postNet:::check_direction("greater"))
    expect_error(postNet:::check_direction(NULL), "cannot be NULL")
    expect_error(postNet:::check_direction(c("greater", "less")), "Please provide only one value")
    expect_error(
        postNet:::check_direction("up"),
        "must be either 'greater' or 'less'"
    )
    expect_null(postNet:::check_category("BP"))
    expect_error(postNet:::check_category(NULL), "cannot be NULL")
    expect_error(postNet:::check_category("Reactome"), "must be a combination")
    expect_null(postNet:::check_analysis_type("rf"))
    expect_error(postNet:::check_analysis_type("glm"), "can only be 'lm' for")
})

test_that("postNetStart validates early user inputs before file work begins", {
    expect_error(postNetStart(source = "load"), "Please provide either an anota2seq object")
    expect_error(
        postNetStart(
            ads = structure(list(), class = "Anota2seqDataSet"),
            geneList = list(groupA = "gene1"),
            source = "load"
        ),
        "Please provide either an anota2seq object or a gene list, not both"
    )
    expect_error(
        postNetStart(
            geneList = list(groupA = "gene1"),
            geneListcolours = c("red", "blue"),
            source = "invalid"
        ),
        "same length as 'geneList'"
    )
    expect_error(
        postNetStart(
            geneList = list(groupA = "gene1"),
            customBg = c("gene2"),
            source = "invalid"
        ),
        "not in 'customBg'"
    )
    expect_error(
        postNetStart(
            geneList = list(groupA = "gene1"),
            source = "invalid"
        ),
        "Invalid source"
    )
})

test_that("additional helper validators cover file, model, feature and motif utilities", {
    expect_true(postNet:::is_by_3(c("ATGAAA", "ATGCCC")))
    expect_false(postNet:::is_by_3(c("ATGA", "ATGCCC")))
    expect_true(postNet:::isKozakContext("adequate1"))
    expect_false(postNet:::isKozakContext("medium"))
    expect_true(postNet:::isUnitOut("position"))
    expect_false(postNet:::isUnitOut("positions"))
    expect_true(postNet:::is_motifs(c("A", "G4")))
    expect_false(postNet:::is_motifs(NULL))

    expect_error(
        postNet:::check_input("createFromSourceFiles", NULL, NULL, "rna.fa.gz", "genomic.gff.gz", NULL, NULL),
        "rna_gbff_file"
    )
    expect_error(
        postNet:::check_input("custom", NULL, NULL, NULL, NULL, NULL, NULL),
        "customFile"
    )
    expect_error(
        postNet:::check_input("createFromFasta", NULL, NULL, NULL, NULL, NULL, NULL),
        "posFile"
    )

    expect_true(postNet:::is_annotType("PTNCDS"))
    expect_false(postNet:::is_annotType("transcript"))
    expect_true(postNet:::is_valid_sourceSeq("create"))
    expect_false(postNet:::is_valid_sourceSeq("remote"))

    tmp_dir <- tempfile("qcdir")
    dir.create(tmp_dir)
    expect_null(postNet:::checkDirectory(tmp_dir))
    expect_error(postNet:::checkDirectory(file.path(tmp_dir, "missing")), "Directory does not exist")

    targetscan_ok <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            Cumulative.weighted.context...score = -0.3,
            Aggregate.PCT = 0.2,
            Gene.Symbol = "gene1",
            Representative.miRNA = "mir1",
            Species.ID = 9606
        ),
        file = targetscan_ok,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    expect_s3_class(postNet:::checkFileColumns(targetscan_ok), "data.frame")

    targetscan_bad <- tempfile(fileext = ".tsv")
    write.table(
        data.frame(
            Cumulative.weighted.context...score = c(-0.3, -0.2),
            Aggregate.PCT = c(0.2, 0.4),
            Gene.Symbol = c("gene1", "gene2"),
            Representative.miRNA = c("mir1", "mir2"),
            Species.ID = c(9606, 10090)
        ),
        file = targetscan_bad,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )
    expect_error(postNet:::checkFileColumns(targetscan_bad), "only the desired species")

    expect_null(postNet:::checkCollection("h"))
    expect_error(postNet:::checkCollection("hallmark"), "valid collections")
    expect_error(postNet:::check_geneList("not-a-list"), "not a list")
    expect_error(postNet:::check_geneList(list()), "list is empty")
    expect_error(postNet:::check_geneList(stats::setNames(list("gene1"), "")), "not a named list")
    expect_null(postNet:::check_size("Count"))
    expect_error(postNet:::check_size(NULL), "must not be NULL")
    expect_error(postNet:::check_size("ratio"), "must not be NULL")
    expect_error(postNet:::check_analysis_type(NULL), "Please provide an input for 'analysis_type'")
    expect_true(postNet:::is_valid_NetModelSel("adjusted"))
    expect_false(postNet:::is_valid_NetModelSel("full"))

    expect_null(postNet:::check_model("finalModel", "lm"))
    expect_null(postNet:::check_model("preModel", "rf"))
    expect_error(postNet:::check_model(NULL, "lm"), "Please provide a valid selection for 'model'")
    expect_error(postNet:::check_model("stepwiseModel", "rf"), "For 'rf', the options are")
    expect_error(postNet:::check_model("borutaModel", "lm"), "valid selection for 'model'")
    expect_null(postNet:::check_features(list(a = 1:2, b = 3:4)))
    expect_error(postNet:::check_features("bad"), "'features' must be a list")
    expect_null(postNet:::check_features(list(a = list(1, 2), b = 3:4)))
    expect_error(postNet:::check_features(list(a = 1:2)), "at least two features")
    expect_null(postNet:::check_lmfeatGroup(c("grp1", "grp2"), 2))
    expect_null(postNet:::check_lmfeatGroup(list("grp1"), 1))
    expect_error(postNet:::check_lmfeatGroup(c("grp1"), 2), "must match the number of 'features'")
    expect_error(postNet:::check_lmfeatGroupColour(c("#AABBCC"), c("grp1")), "has no names")
    expect_null(postNet:::check_lmfeatGroupColour(c(grp1 = "#AABBCC", grp2 = "#112233"), c("grp1", "grp2")))
    expect_error(postNet:::check_lmfeatGroupColour(c(grp1 = "#AABBCC"), c("grp1", "grp2")), "do not exactly match")
    expect_error(postNet:::check_lmfeatGroupColour(c(grp1 = "red"), c("grp1")), "hex color")
    expect_true(postNet:::check_shiftUnit("p10"))
    expect_false(postNet:::check_shiftUnit("P10"))

    feats <- data.frame(a = 1:2, b = 2:1, row.names = c("g1", "g2"))
    expect_true(postNet:::check_featSel(c("a", "b"), feats))
    expect_false(postNet:::check_featSel("a", feats))
    expect_true(postNet:::check_featCol("a", feats))
    expect_false(postNet:::check_featCol("missing", feats))
    expect_true(postNet:::check_predFeat(data.frame(a = c(1, 2), b = c(2, 3), row.names = c("g1", "g2"))))
    expect_false(postNet:::check_predFeat(data.frame(a = c("x", "y"), row.names = c("g1", "g2"))))
    expect_true(postNet:::check_predFeat(structure(data.frame(a = c(1, 2)), row.names = NULL)))
    expect_false(postNet:::check_predFeat(data.frame(a = c(1, NA), row.names = c("g1", "g2"))))
})
