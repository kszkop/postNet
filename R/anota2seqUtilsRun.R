anota2seqUtilsRun <- function(source, 
                              species, 
                              version=NULL,
                              
                              region, 
                              ads, 
                              contrast, 
                              regulation, 
                              comparisons, 
                              selection, 
                              plotType, 
                              pdfName){
  
  #
  if(is.null(version)){
    version <- checkAvailableVersions(species=species)
  }
  #
  annot <- retrieveFormatData(source = source, 
                              version = version,
                              species = species)
  #
  for(i in region){
    
    lengthAnalysis(ads=ads, 
                  annot=annot, 
                  regulation=regulation, 
                  contrast=contrast, 
                  comparisons=comparisons, 
                  region=i, 
                  selection=selection, 
                  plotType=plotType, 
                  pdfName=pdfName)
    
    contentAnalysis(ads = ads, 
                    annot = annot, 
                    regulation =  regulation,
                    contrast = contrast, 
                    region = i, 
                    selection = selection,
                    comparisons = comparisons, 
                    content = content,
                    plotType = plotType, 
                    pdfName = pdfName)
    for(j in 1:length(regulation)){
      motifAnalysis(ads = ads, 
                    annot = annot, 
                    regulation = regulation[i], 
                    contrast = contrast[i],
                    region = i, 
                    selection = selection)
    }
    
    ##combine all motifs for given region
    motifsSel <- c(motifs_utr5_transUp$consensus,motifs_utr5_transDown$consensus)
    motifsOut <- lapply(motifsSel, function(x) contentMotifs(ads = ads, annot = annot, regulation =  c("translationUp","translationDown"),contrast = c(1,1), motif=x,region = 'UTR5', len=1, selection ='random',comparisons = list(c(1,2))))
    names(motifsOut) <- paste('UTR5', motifsSel, sep='_')
    
    if(i != 'UTR3'){
        contentMotifs(ads = ads, 
                      annot = annot, 
                      regulation =  c("translationUp","translationDown"), 
                      contrast = c(1,1), 
                      motif='G4',
                      region = 'UTR5', 
                      len=1, 
                      selection ='random', 
                      comparisons = list(c(1,2)))
    
    }
    
    foldEnergyAnalysis(ads = ads, 
                       annot = annot, 
                       regulation = regulation, 
                       contrast =  ccontrast, 
                       region= i,
                       selection = selection, 
                       comparisons = comparisons, 
                       source = source, 
                       version = version,
                       species = species, 
                       resid = resid, 
                       onlyRun = onlyRun, 
                       plotType = plotType, 
                       pdfName = NULL)
    
    
    
    
    
  }
  
  uorf_analysis(ads = ads, 
                annot = annot, 
                regulation =  regulation, 
                contrast = contrast, 
                selection ='random', 
                comparisons = list(c(0,1),c(0,2),c(1,2)), 
                onlyUTR5=FALSE, 
                startCodon='ATG', 
                KozakContext='strong')
  
  
  
}