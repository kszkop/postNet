slopeFilt <- function(ads,
                      regulationGen,
                      contrastSel,
                      minSlope, 
                      maxSlope ){
    #
    tmpAds <- anota2seq::anota2seqGetOutput(ads, 
                                 analysis = regulationGen,
                                 output = "full",
                                 selContrast = contrastSel,
                                 getRVM = TRUE)
    
    #Filter slopes
    tmpAds_slopeFilt <- tmpAds[which(tmpAds[,1] < minSlope | tmpAds[,1] > maxSlope),]
    
    #vector of genes to out
    genesOut <- as.character(row.names(tmpAds_slopeFilt))
    #
    return(genesOut)
}
