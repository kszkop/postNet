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