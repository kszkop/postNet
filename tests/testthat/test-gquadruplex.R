test_that("find_gquadruplexes returns positions and counts for canonical motifs", {
    hits <- postNet:::find_gquadruplexes("GGGTGGGTGGGTGGG", min_score = 26, return = "position")
    expect_s3_class(hits, "data.frame")
    expect_equal(hits$start, 1)
    expect_equal(hits$end, 15)

    count <- postNet:::find_gquadruplexes("GGGTGGGTGGGTGGG", min_score = 47, return = "number")
    expect_equal(count, 1L)
})

test_that("find_gquadruplexes returns empty outputs when no G4 is found", {
    hits <- postNet:::find_gquadruplexes("ACTTACTTACTT", return = "position")
    expect_s3_class(hits, "data.frame")
    expect_equal(nrow(hits), 0)

    count <- postNet:::find_gquadruplexes("ACTTACTTACTT", return = "number")
    expect_equal(count, 0L)
})

test_that("calc_g4 maps helper output to number and position forms", {
    expect_equal(postNet:::calc_g4("GGGTGGGTGGGTGGG", min_score = 26, unit = "number"), 1)

    pos <- postNet:::calc_g4("GGGTGGGTGGGTGGG", min_score = 26, unit = "position")
    expect_equal(pos$start, 1)
    expect_equal(pos$end, 15)

    expect_equal(postNet:::calc_g4("ACTTACTTACTT", min_score = 26, unit = "number"), 0)
    expect_true(is.na(postNet:::calc_g4("ACTTACTTACTT", min_score = 26, unit = "position")))
})
