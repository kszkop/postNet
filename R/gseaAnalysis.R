gseaAnalysis <- function(
        ptn,
        genesSlopeFiltOut = NULL,
        collection = NULL,
        subcollection = NULL,
        subsetNames = NULL,
        geneSet = NULL,
        maxSize = 500,
        minSize = 10,
        name = NULL
) {
    check_ptn(ptn)
    
    if (is.null(geneSet) && is.null(collection)) {
        stop("Please provide an input for either 'geneSet' or 'collection'.")
    }
    
    if (!check_number(maxSize) || !check_number(minSize)) {
        stop("The inputs for 'maxSize' and 'minSize' must be integers.")
    }
    
    if (minSize <= 0 || maxSize <= 0) {
        stop("The inputs for 'maxSize' and 'minSize' must be positive.")
    }
    
    if (maxSize <= minSize) {
        stop("'maxSize' must be greater than 'minSize'.")
    }
    
    effTmp <- ptn_effect(ptn)
    
    if (!is.null(genesSlopeFiltOut)) {
        effIn <- effTmp[!names(effTmp) %in% genesSlopeFiltOut]
    } else {
        effIn <- effTmp
    }
    
    rankIn <- effIn[order(effIn, decreasing = TRUE)]
    
    if (length(rankIn) == 0) {
        stop("No ranked genes are available for GSEA after filtering.")
    }
    
    if (is.null(geneSet)) {
        species <- ptn_species(ptn)
        
        if (!species %in% c("human", "mouse")) {
            stop(
                "This option is currently only available for human or mouse."
            )
        }
        
        checkCollection(collection)
        
        eh <- ExperimentHub::ExperimentHub()
        AnnotationHub::query(eh, "msigdb")
        
        versionTmp <- as.character(
            sort(
                as.numeric(msigdb::getMsigdbVersions()),
                decreasing = TRUE
            )[1]
        )
        
        msigdbOut <- msigdb::getMsigdb(
            org = ifelse(species == "human", "hs", "mm"),
            id = "SYM",
            version = versionTmp
        )
        
        msigdbOut <- msigdb::appendKEGG(msigdbOut, version = versionTmp)
        
        collectionTmp <- msigdb::subsetCollection(
            msigdbOut,
            collection = collection,
            subcollection = subcollection
        )
        
        if (!is.null(subsetNames)) {
            collectionTmp <- collectionTmp[
                names(collectionTmp) %in% subsetNames
            ]
        }
        
        geneSet_ids <- GSEABase::geneIds(collectionTmp)
    } else {
        check_geneList(geneSet)
        geneSet_ids <- geneSet
    }
    
    if (length(geneSet_ids) == 0) {
        stop("No gene sets are available for analysis.")
    }
    
    resOut <- fgsea::fgsea(
        pathways = geneSet_ids,
        stats = rankIn,
        minSize = minSize,
        maxSize = maxSize
    )
    
    nameTmp <- ifelse(
        !is.null(name),
        paste(name, "gseaAnalysis", sep = "_"),
        "gseaAnalysis"
    )
    
    if (nrow(resOut) == 0) {
        gseaOut <- data.frame(
            Term = character(),
            ES = numeric(),
            NES = numeric(),
            log2err = numeric(),
            Count = integer(),
            Size = integer(),
            pvalue = numeric(),
            adjusted_pvalue = numeric(),
            Genes = character(),
            stringsAsFactors = FALSE
        )
        
        data.table::fwrite(
            gseaOut,
            file = paste(nameTmp, ".txt", sep = ""),
            sep = "\t"
        )
        
        ptn@analysis@GSEA <- gseaOut
        
        return(ptn)
    }
    
    resOut$Count <- vapply(resOut$leadingEdge, length, integer(1))
    
    colnames(resOut) <- c(
        "Term",
        "pvalue",
        "adjusted_pvalue",
        "log2err",
        "ES",
        "NES",
        "Size",
        "Genes",
        "Count"
    )
    
    resOut <- resOut[, c(1, 5, 6, 4, 9, 7, 2, 3, 8)]
    
    gseaOut <- resOut[order(resOut$adjusted_pvalue), ]
    
    gseaOut$Genes <- vapply(
        gseaOut$Genes,
        function(x) {
            paste(x, collapse = ":")
        },
        character(1)
    )
    
    data.table::fwrite(
        gseaOut,
        file = paste(nameTmp, ".txt", sep = ""),
        sep = "\t"
    )
    
    ptn@analysis@GSEA <- gseaOut
    
    return(ptn)
}
