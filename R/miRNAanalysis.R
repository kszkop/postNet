miRNAanalysis <- function(
        ptn,
        miRNATargetScanFile,
        genesSlopeFiltOut = NULL,
        contextScore = -0.2,
        Pct = 0,
        maxSize = 500,
        minSize = 10
) {
    check_ptn(ptn)
    
    if (!check_number(maxSize) || !check_number(minSize)) {
        stop("Please provide numeric values for 'maxSize' and 'minSize'.")
    }
    
    if (minSize <= 0 || maxSize <= 0) {
        stop("The size parameters must be positive integers.")
    }
    
    if (maxSize <= minSize) {
        stop("maxSize must be greater than minSize.")
    }
    
    if (is.null(contextScore) && is.null(Pct)) {
        stop(
            "Please provide values for at least one of 'contextScore' ",
            "or 'Pct'."
        )
    }
    
    if (!is.null(Pct)) {
        if (!check_number(Pct)) {
            stop("Please provide a numeric value for 'Pct'.")
        }
        
        if (Pct < 0 || Pct > 1) {
            stop("Please provide a value between 0 and 1 for 'Pct'.")
        }
    }
    
    if (!is.null(contextScore)) {
        if (!check_number(contextScore)) {
            stop("Please provide a numeric value for 'contextScore'.")
        }
        
        if (contextScore > 0) {
            stop("Please provide a negative value for 'contextScore'.")
        }
    }
    
    if (!is.null(Pct) && !is.null(contextScore)) {
        if (Pct == 0 && contextScore == 0) {
            message(
                "Consider providing filtering thresholds for ",
                "'contextScore' or 'Pct'."
            )
        }
    }
    
    miRNATargetScan <- checkFileColumns(miRNATargetScanFile)
    
    if (
        length(intersect(
            unique(miRNATargetScan$Gene.Symbol),
            ptn_background(ptn)
        )) < 2
    ) {
        stop(
            "No overlap detected between genes in the dataset and in ",
            "the miRNATargetScanFile."
        )
    }
    
    effTmp <- ptn_effect(ptn)
    
    if (!is.null(genesSlopeFiltOut)) {
        effIn <- effTmp[!names(effTmp) %in% genesSlopeFiltOut]
    } else {
        effIn <- effTmp
    }
    
    rankIn <- effIn[order(effIn, decreasing = TRUE)]
    
    if (!is.null(contextScore)) {
        miRNATargetScan <- miRNATargetScan[
            miRNATargetScan$Cumulative.weighted.context...score <
                contextScore,
        ]
    }
    
    if (!is.null(Pct)) {
        miRNATargetScan <- miRNATargetScan[
            miRNATargetScan$Aggregate.PCT > Pct,
        ]
    }
    
    miRNA_names <- unique(miRNATargetScan$Representative.miRNA)
    
    if (length(miRNA_names) == 0) {
        stop(
            "No miRNA meet the selected thresholds. Try adjusting ",
            "'contextScore' and/or 'Pct' for less stringent filtering."
        )
    }
    
    miRNAList <- vector("list", length(miRNA_names))
    names(miRNAList) <- miRNA_names
    
    for (i in seq_along(miRNAList)) {
        miRNAList[[i]] <- miRNATargetScan$Gene.Symbol[
            miRNATargetScan$Representative.miRNA %in% names(miRNAList)[i]
        ]
    }
    
    miRNA_enrichment_targetScan <- gage::gage(
        exprs = rankIn,
        gsets = miRNAList,
        set.size = c(minSize, maxSize),
        rank.test = TRUE,
        use.fold = FALSE
    )
    
    miRNA <- new(
        "postNetmiRNA",
        miRNA_analysis = miRNA_enrichment_targetScan,
        miRNA_to_gene = miRNAList
    )
    
    ptn@analysis@miRNA <- miRNA
    
    return(ptn)
}
