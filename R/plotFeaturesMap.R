plotFeaturesMap <- function(ptn,
                            regOnly = TRUE,
                            comparisons = NULL,
                            featSel,
                            remBinary = TRUE,
                            featCol = NULL,
                            scaled = FALSE,
                            centered = TRUE,
                            remExtreme = NULL,
                            pdfName = NULL) {
  check_ptn(ptn)

  if (!is.null(ptn_features(ptn))) {
    featuresIn <- ptn_features(ptn)
  } else {
    stop("Features analysis is NULL. Please run the featureIntegration() analysis first.")
  }
  if (!check_featSel(featSel, features = featuresIn)) {
    stop("Please provide a character vector of features of interest to be used for UMAP embeddings. We recommed starting with those that were selected in the \ final featureIntegration() models. These can be obtained by running the ptn_selectedFeatures() function. See help manuals and vignette for details.")
  }
  if (!is.null(featCol)) {
    if (!check_featCol(featCol, features = featuresIn)) {
      stop("Please provide a character vector defining the features that will be overlaid on UMAPs. All features must be present in the 'features' stored in \ the ptn object. To check which features have been stored run: colnames(ptn_features(ptn))")
    }
  }
  featuresSel <- featuresIn[, colnames(featuresIn) %in% featSel]

  if (!check_logical(regOnly)) {
    stop("The input for 'regOnly' must be logical: TRUE or FALSE")
  }
  if (!check_logical(scaled)) {
    stop("The input for 'scaled' must be logical: TRUE or FALSE")
  }
  if (!check_logical(centered)) {
    stop("The input for 'centered' must be logical: TRUE or FALSE")
  }
  if (!check_logical(remBinary)) {
    stop("The input for 'remBinary' must be logical: TRUE or FALSE")
  }
  if (isTRUE(regOnly)) {
    if (!check_comparisons(comparisons)) {
      stop("The input for 'comparisons' must be a list of numeric vectors of paired comparisons. For example: list(c(0,2), c(0,1)). 0 always \ denotes the background gene set.")
    }
    if (length(which(unique(unlist(comparisons)) == 0)) > 0 && is.null(ptn_background(ptn))) {
      stop("0 always denotes the background, but no background has been provided.")
    }
    if (length(comparisons) != 1) {
      stop("Although potentially possible, please run each comparison separately to fit your feature selection.")
    }
  }
  if (!is.null(remExtreme)) {
    if (!check_number(remExtreme)) {
      stop("The input for 'remExtreme' should be a number between 0 and 1 indicating a percentile threshold for extreme values to be removed \ for determining the colouring scale.")
    } else {
      if (remExtreme <= 0 | remExtreme >= 1) {
        stop("The input for 'remExtreme' can only be a single value between 0 and 1. It represents a percentile threshold in decimal format for extreme values \ from both sides of the feature distribution to be removed. Values are removed only for determining the colouring scale.")
      }
    }
  }
  #
  if (length(which(unique(unlist(comparisons)) == 0)) > 0 && is.null(ptn_background(ptn))) {
    stop("0 always denotes the background, but no background has been provided.")
  }
  #
  if (isTRUE(regOnly)) {
    resOut <- resQuant(qvec = ptn_effect(ptn), ptn = ptn)
    # for (i in 1:length(comparisons)) {
    if (names(resOut)[1] == "background") {
      compTmp <- comparisons[[1]] + 1
    } else {
      compTmp <- comparisons[[1]]
    }
    listSel <- c(names(resOut[[compTmp[1]]]), names(resOut[[compTmp[2]]]))
    featuresSel <- featuresSel[row.names(featuresSel) %in% listSel, ]
  }
  if (isTRUE(remBinary)) {
    binaryCols <- sapply(featuresSel, is_binary)
    featuresCluster <- featuresSel[, !binaryCols, drop = FALSE]
  } else {
    featuresCluster <- featuresSel
  }
  if (isTRUE(scaled)) {
    # save <- featuresCluster
    featuresCluster <- scale(na.omit(featuresCluster), center = centered, scale = scaled)
  } else {
    featuresCluster <- na.omit(featuresCluster)
  }

  fmapRes <- umap::umap(featuresCluster, n_components = 2)
  fmapRes <- fmapResOut <- as.data.frame(fmapRes$layout)
  colnames(fmapRes) <- colnames(fmapResOut) <- c("UMAP1", "UMAP2")
  fmapRes$Gene <- rownames(featuresCluster)

  # Plot selected features

  effTmp <- ptn_effect(ptn)
  eff <- effTmp[match(row.names(fmapRes), names(effTmp))]
  effect_fmap <- plot_fmap(fmapRes, colVec = eff, remExtreme = remExtreme, name = "Effect")

  if (is.null(featCol)) {
    featCol <- featSel
  }

  featuresIn <- featuresIn[match(row.names(fmapRes), row.names(featuresIn)), ]
  for (feat in featCol) {
    featTmp <- featuresIn[, feat]
    feature_fmap <- plot_fmap(fmapRes, colVec = featTmp, remExtreme = remExtreme, name = feat)

    pdf(ifelse(is.null(pdfName), paste(feat, "_featureUMAP.pdf", sep = ""), paste(pdfName, feat, "featureUMAP.pdf", sep = "_")), width = 16, height = 8, useDingbats = FALSE)
    par(mar = c(5, 5, 5, 5), bty = "l", font = 2, font.axis = 2, font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.9)

    gridExtra::grid.arrange(
      effect_fmap$legend[[1]],
      effect_fmap$mainPlot,
      feature_fmap$mainPlot,
      feature_fmap$legend[[1]],
      ncol = 4,
      widths = c(1, 4, 4, 1)
    )

    dev.off()
  }
  ptn@analysis@featureIntegration[["featuresMap"]] <- fmapResOut
  return(ptn)
}
