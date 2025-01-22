slopeFilt <- function(ads,
                      regulationGen,
                      contrastSel,
                      minSlope, 
                      maxSlope ){
  #
  if (!check_ads(ads)) {
    stop("ads is not a valid 'Anota2seqDataSet' object.")
  }
  if (!regulationGen %in% c("translation","buffering")) {
    stop("For filtering slopes, 'regulationGen' should be a character vector chosen from 'translation' or 'buffering' ")
  }
  if(!is.numeric(contrastSel) | !contrastSel %in% seq(1,ncol(ads@contrasts),1)){
    stop("'contrastSel' should be a numeric for the desired contrast in anota2seq object")
  }
  checkSlopes(minSlope, maxSlope)
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

