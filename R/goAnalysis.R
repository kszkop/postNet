goAnalysis <- function(ptn,
                       genesSlopeFiltOut = NULL,
                       category,
                       maxSize = 500,
                       minSize = 10,
                       counts = 10,
                       FDR = 0.15,
                       name = NULL) {
  #
  check_ptn(ptn)
  species <- ptn_species(ptn)
  if (!species %in% c("human", "mouse")) {
    stop("GO term analysis functionality is currently only available for human and mouse.")
  }
  #
  check_category(category)
  #
  if (!check_number(maxSize) | !check_number(minSize) | !check_number(counts) | !check_number(FDR)) {
    stop("The inputs for 'maxSize', 'minSize', 'counts', and 'FDR' must be numeric.")
  }
  if (minSize <= 0 | maxSize <= 0 | counts <= 0 | FDR < 0) {
    stop("The inputs for 'maxSize', 'minSize', 'counts', and 'FDR' must be positive.")
  }
  if (maxSize <= minSize) {
    stop("'maxSize' must be greater than 'minSize'.")
  }
  #
  GOout <- list()

  res <- ptn_geneList(ptn)
  bg <- unlist(ptn_background(ptn))
  if (length(setdiff(bg, unique(unlist(res)))) == 0) {
    warning("All genes in the background gene set are regulated. Please check that you are using an appropriate background set.")
  }

  # filter for slopes if indicated
  if (!is.null(genesSlopeFiltOut)) {
    bg <- bg[!bg %in% genesSlopeFiltOut]
    res <- lapply(res, function(x) x[!x %in% genesSlopeFiltOut])
  }

  # Convert to Entrez ID
  bg_entrezID <- convertSymbolToEntrezID(geneList = bg, species = species)
  res_entrezID <- lapply(res, function(x) convertSymbolToEntrezID(geneList = x, species = species))

  GoLists <- res_entrezID[!sapply(res_entrezID, is.null)]

  #
  GOout <- new("postNetGO",
    BP = NULL,
    CC = NULL,
    MF = NULL,
    KEGG = NULL
  )
  #
  for (sel in category) {
    resOut <- list()
    if (sel == "KEGG") {
      for (i in 1:length(GoLists)) {
        resTmp <- clusterProfiler::enrichKEGG(
          gene = GoLists[[i]],
          universe = bg_entrezID,
          organism = ifelse(species == "human", "hsa", ifelse(species == "mouse", "mmu", "ups")),
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          minGSSize = minSize,
          maxGSSize = maxSize
        )
        #
        resOut[[names(GoLists)[i]]] <- resTmp
      }
    } else if (sel == "BP" | sel == "MF" | sel == "CC") {
      for (i in 1:length(GoLists)) {
        resTmp <- clusterProfiler::enrichGO(
          gene = GoLists[[i]],
          universe = bg_entrezID,
          OrgDb = ifelse(species == "human", "org.Hs.eg.db", ifelse(species == "mouse", "org.Mm.eg.db", "ups")),
          ont = sel,
          pAdjustMethod = "BH",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          minGSSize = minSize,
          maxGSSize = maxSize,
          readable = TRUE,
          pool = FALSE
        )
        resOut[[names(GoLists)[i]]] <- resTmp
      }
    }
    # remove counts below 10 and recalculate adjp and reformat
    for (i in 1:length(resOut)) {
      tabTmp <- resOut[[i]]@result
      #
      tabTmp <- tabTmp[tabTmp$Count >= counts, ]
      if (nrow(tabTmp) > 0) {
        tabTmp$p.adjust <- stats::p.adjust(tabTmp$pvalue, method = "BH")
      }
      tabTmp <- tabTmp[tabTmp$p.adjust < FDR, ]

      if (nrow(tabTmp) > 0) {
        geneIDs_temp <- tabTmp$geneID

        checkID <- check_id_type(seqinr::c2s(strsplit(geneIDs_temp[[1]], "/")[[1]][1:5]))
        if (checkID == "entrezID") {
          tabTmp$geneID <- sapply(geneIDs_temp, function(x) paste(sort(convertEntrezIDToSymbol(unlist(strsplit(x, "/")), species = species)), collapse = ":"), USE.NAMES = F)
        } else {
          tabTmp$geneID <- sapply(geneIDs_temp, function(x) paste(sort(unlist(strsplit(x, "/"))), collapse = ":"), USE.NAMES = F)
        }
      } else {
        message(paste("No significant results for", sel, "in", names(resOut)[i], sep = " "))
      }
      tabTmp$Size <- as.numeric(sub("\\/.*", "", tabTmp$BgRatio))

      if (sel == "KEGG") {
        tabTmp <- tabTmp[, c("ID", "Description", "category", "subcategory", "Count", "Size", "pvalue", "p.adjust", "geneID")]
      } else {
        tabTmp <- tabTmp[, c("ID", "Description", "Count", "Size", "pvalue", "p.adjust", "geneID")]
      }
      resOut[[i]]@result <- tabTmp
    }
    resWrite <- lapply(resOut, function(x) x@result)
    #
    nameOut <- ifelse(is.null(name), paste("GO_", sel, ".xlsx", sep = ""), paste(name, "_GO_", sel, ".xlsx", sep = ""))
    WriteXLS::WriteXLS(resWrite, SheetNames = names(resWrite), ExcelFileName = nameOut, row.names = FALSE)

    slot(GOout, sel) <- resOut
  }
  ptn@analysis@GO <- GOout
  #
  return(ptn)
}
