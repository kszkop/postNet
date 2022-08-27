
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
  if(is.null(version){
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
    
    
  }
  
  
  
  
  
}