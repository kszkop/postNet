slopeFilt <- function(ads,
                      regulation, #'translation' or 'buffering'
                      contrast,
                      minSlope, 
                      maxSlope
){
    #extract results from anota2seq
    tmpAds <- anota2seq::anota2seqGetOutput(ads, 
                                 analysis = regulation,
                                 output = "full",
                                 selContrast = contrast,
                                 getRVM = TRUE)
    
    #Filter slopes
    tmpAds_slopeFilt <- tmpAds[which(tmpAds[,1] < minSlope | tmpAds[,1] > maxSlope),]
    
    #vector of genes to out
    genesOut <- as.character(row.names(tmpAds_slopeFilt))
    
    #
    return(genesOut)
}
