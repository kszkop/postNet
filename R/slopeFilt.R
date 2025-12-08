slopeFilt <- function(ads,
                      regulationGen,
                      contrastSel,
                      minSlope = NULL,
                      maxSlope = NULL) {
  #
  check_ads(ads)
  if (!regulationGen %in% c("translation", "buffering")) {
    stop("For filtering slopes, 'regulationGen' should be either 'translation' or 'buffering'.")
  }
  if (!is.numeric(contrastSel) | !contrastSel %in% seq(1, ncol(ads@contrasts), 1)) {
    stop("THe input for 'contrastSel' should be a number corresponding to the desired contrast in anota2seq object.")
  }

  # Set default values for minSlope and maxSlope based on regulationGen if they are NULL
  if (is.null(minSlope) || is.null(maxSlope)) {
    if (regulationGen == "translation") {
      minSlope <- ifelse(is.null(minSlope), -1, minSlope)
      maxSlope <- ifelse(is.null(maxSlope), 2, maxSlope)
    } else if (regulationGen == "buffering") {
      minSlope <- ifelse(is.null(minSlope), -2, minSlope)
      maxSlope <- ifelse(is.null(maxSlope), 1, maxSlope)
    }
  }

  checkSlopes(minSlope, maxSlope)
  #
  tmpAds <- anota2seq::anota2seqGetOutput(ads,
    analysis = regulationGen,
    output = "full",
    selContrast = contrastSel,
    getRVM = TRUE
  )

  # Filter slopes
  tmpAds_slopeFilt <- tmpAds[which(tmpAds[, 1] < minSlope | tmpAds[, 1] > maxSlope), ]

  # vector of genes to out
  genesOut <- as.character(row.names(tmpAds_slopeFilt))
  #
  return(genesOut)
}
