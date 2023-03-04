anota2seqUtilsRun <- function(source='load', 
                              version=NULL,
                              species, 
                              customFile=NULL,
                              fastaFile=NULL, 
                              posFile=NULL, 
                              rna_gbff_file=NULL, 
                              rna_fa_file=NULL, 
                              genomic_gff_file=NULL,
                              
                              adjObj=NULL,
                              region_adj=NULL,
                              excl=FALSE,
                              keepAll=FALSE,
                              
                              regionSel,
                              
                              ads = NULL,
                              regulation = NULL,
                              contrast = NULL,
                              geneList = NULL,
                              geneListcolours = NULL,
                              customBg = NULL,
                              selection = "random",
                              comparisons = NULL,
                              plotOut = TRUE,
                              plotType = "boxplot",
                              pdfName = NULL,
   
                              contentIn,
                              subregion=NULL, 
                              subregionSel=NULL, 
                              
                              startCodon='ATG',
                              KozakContext='strong',
                              onlyUTR5=FALSE,
                              
                              runMotifs=TRUE,
                              stremeThreshold = 0.05,
                              minwidth = 6,
                              memePath = NULL,
                              
                              seqType='dna', 
                              dist=1,
                              min_score=47, 
                              resid=FALSE,
                              
                              sourceFE='load',
                              fromFasta=FALSE,
                              customFile=NULL, 
                              onlyRun=FALSE,
                              residFE=FALSE,
                              
                              addMotifs=NULL,
                              addSeqType=NULL,
                    
                              annotType = "ccds",
                              sourceCod = "load",
                              customFileCod = NULL,
                              analysis="codon",
                              codonN = 1,
                              type = "sequence",
                              dpn_path = NULL,
                              cds_filt = TRUE,
                              pAdj = 0.01,
                              plotHeatmap = TRUE,
                              thresX1 = 0.3,
                              thresY1 = 0.3,
                              thresX2 = 0.3,
                              thresY2 = 0.3,
                              
                              unit = "count",
                              
                              addSign = NULL,
                              
                              regulationGen = NULL,
                              minSlope = NULL,
                              maxSlope = NULL,
                              
                              miRNAmax=25,
                              miRNATargetScanFile=NULL,
                              
                              analysis_type,
                              rmCat = NULL,
                              regOnly = TRUE,
                              allFeat = TRUE,
                              covarFilt = 20,
                              NetModelSel = "Omnibus",
                              effectMeasure = NULL,
                              outDir = NULL
                              ){
  #####Prepare annotation
  if(is.null(version)){
    version <- checkAvailableVersions(species=species)
    #extract the latest
    versionInd <- sub("^[^.]*.","", version)
    versionInd <- sort(versionInd,decreasing = T)[1]
    version <- version[grep(versionInd, version)]
  }
  #
  annot <- retrieveFormatData(source=source, 
                              species=species,
                              version=version,
                              customFile=customFile, 
                              fastaFile=fastaFile, 
                              posFile=posFile,
                              rna_gbff_file=rna_gbff_file, 
                              rna_fa_file=rna_fa_file, 
                              genomic_gff_file=genomic_gff_file)
  
  if(!is.null(adjObj)){
    annot <- adjustSeq(annot=annot, 
                       adjObj = adjObj,
                       region_adj = region_adj,
                       excl = excl,
                       keepAll = keepAll)
  }   
  #####Run analysis
  resultsOut <- list()
  #
  for(i in 1:regionSel){
    #Select current region to run
    region <- regionSel[i]
    #
    resultsOut[[paste('len',region,sep='_')]] <- lengthAnalysis(annot=annot, 
                                                                ads=ads,
                                                                regulation=regulation, 
                                                                contrast=contrast, 
                                                                geneList=geneList,
                                                                geneListcolours=geneListcolours,
                                                                customBg=customBg,
                                                                selection=selection,
                                                                region=region, 
                                                                comparisons=comparisons,
                                                                plotOut=plotOut,
                                                                plotType=plotType, 
                                                                pdfName=pdfName)
    #
    contentTmp <- contentAnalysis(annot=annot,
                                  contentIn=contentIn,
                                  ads = ads,
                                  regulation = regulation,
                                  contrast = contrast,
                                  geneList = geneList,
                                  geneListcolours = geneListcolours,
                                  customBg = customBg,
                                  selection = selection,
                                  region = region,
                                  comparisons = comparisons,
                                  plotOut = plotOut,
                                  plotType = plotType,
                                  pdfName = pdfName)
    resultsOut <- append(resultsOut, contentTmp)
    #
    if(region=='UTR5'){
      resultsOut[[paste('uORFs',startCodon,KozakContext,sep='_')]] <- uorf_analysis(annot = annot,
                                                                                    startCodon=startCodon,
                                                                                    KozakContext=KozakContext,
                                                                                    onlyUTR5=onlyUTR5,
                                                                                    unitOut="number",                                                                                    ads = ads, 
                                                                                    regulation =  regulation, 
                                                                                    contrast = contrast,
                                                                                    geneList = geneList,
                                                                                    geneListcolours = geneListcolours,
                                                                                    customBg = customBg,
                                                                                    selection = selection,
                                                                                    comparisons = comparisons,
                                                                                    plotOut = plotOut,
                                                                                    pdfName = pdfName)
    }
    #
    if(isTRUE(runMotifs)){
      tmpOut <- motifAnalysis(annot = annot,
                              stremeThreshold = stremeThreshold,
                              minwidth=minwidth,
                              memePath=memePath,
                              ads = ads, 
                              regulation=regulation,
                              contrast = contrast,
                              geneList=geneList,
                              customBg=customBg,
                              selection = selection,
                              region = region,
                              subregion=subregion, 
                              subregionSel= subregionSel)
      motifsIn <- tmpOut[[1]]
      #
      if(length(motifsIn)>0){
        #
        motifsOut <- contentMotifs(annot=annot,
                                   motifsIn=motifsIn,
                                   seqType = seqType,
                                   dist = dist,
                                   min_score = min_score,
                                   resid = resid,
                                   ads = ads,
                                   regulation = regulation,
                                   contrast = contrast,
                                   geneList = geneList,
                                   geneListcolours = geneListcolours,
                                   customBg = customBg,
                                   selection=selection,
                                   region=region,
                                   subregion = subregion,
                                   subregionSel = subregionSel,
                                   comparisons = comparisons,
                                   pdfName = pdfName,
                                   plotOut = plotOut)
        
        resultsOut <- append(resultsOut, motifsOut)
      } else {
        print('No significant de-novo motifs found')
      }
    }
    if(i != 'UTR3'){
      resultsOut[[paste(region,'G4',sep='_')]] <- contentMotifs(annot = annot,
                                                                motifsIn = "G4",
                                                                seqType = seqType,
                                                                dist = dist,
                                                                min_score = min_score,
                                                                resid = resid,
                                                                ads = ads,
                                                                regulation = regulation,
                                                                contrast = contrast,
                                                                geneList = geneList,
                                                                geneListcolours = geneListcolours,
                                                                customBg = customBg,
                                                                selection=selection,
                                                                region=region,
                                                                subregion = subregion,
                                                                subregionSel = subregionSel,
                                                                comparisons = comparisons,
                                                                pdfName = pdfName,
                                                                plotOut = plotOut)[[1]]
    }
    #
    if(!is.null(addMotifs)){
      motifsOutadd <- contentMotifs(annot = annot,
                                    motifsIn = addMotifs,
                                    seqType = addSeqType,
                                    dist = dist,
                                    min_score = min_score,
                                    resid = resid,
                                    ads = ads,
                                    regulation = regulation,
                                    contrast = contrast,
                                    geneList = geneList,
                                    geneListcolours = geneListcolours,
                                    customBg = customBg,
                                    selection=selection,
                                    region=region,
                                    subregion = subregion,
                                    subregionSel = subregionSel,
                                    comparisons = comparisons,
                                    pdfName = pdfName,
                                    plotOut = plotOut)
      
      resultsOut <- append(resultsOut, motifsOutadd)
    }
    #
    resultsOut[[paste('fe',region,sep='_')]] <- foldingEnergyAnalysis(annot = annot, 
                                                                      sourceFE = sourceFE,
                                                                      version=version,
                                                                      species=species,
                                                                      fromFasta = fromFasta,
                                                                      customFile = customFile,
                                                                      onlyRun = onlyRun, 
                                                                      residFE = residFE,
                                                                      ads = ads,
                                                                      regulation = regulation,
                                                                      contrast = contrast,
                                                                      geneList = geneList,
                                                                      geneListcolours = geneListcolours,
                                                                      customBg = customBg,
                                                                      selection=selection,
                                                                      region=region,
                                                                      comparisons = comparisons,
                                                                      plotOut = plotOut,
                                                                      plotType = plotType,
                                                                      pdfName = pdfName)
    #
    if(region=='CDS'){
      codonOutTmp <- codonUsage(annot=annot,
                                annotType = annotType, 
                                sourceCod = sourceCod,
                                customFileCod = customFile,
                                species = species,
                                analysis=analysis,
                                codonN = codonN,
                                type = type, 
                                dpn_path = dpn_path,
                                cds_filt = cds_filt,
                                pAdj =  pAdj,
                                plotHeatmap = plotHeatmap,
                                thresX1 = thresX1,
                                thresY1 = thresY1,
                                thresX2 = thresX2,
                                thresY2 = thresY2,
                                ads = ads,
                                regulation = regulation,
                                contrast = contrast,
                                geneList = geneList,
                                geneListcolours = geneListcolours,
                                customBg = customBg,
                                selection=selection, 
                                subregion =subregion, 
                                subregionSel=subregionSel,
                                comparisons=comparisons,
                                pdfName = pdfName)
      
      #
      codonIn <- codonOutTmp[['codonAll']]
      for(j in 1:length(comparisons)){
        #
        featSelUp <- codonOutTmp[[j]][["Up"]]
        featSelDown <- codonOutTmp[[j]][["Down"]]
        #
        if(length(featSelUp) > 0){
          selCodonOutUp <- codonCalc(codonIn = codonIn,
                                     featsel = featSelUp,
                                     analysis = analysis,
                                     unit = unit,
                                     ads = ads,
                                     regulation = regulation,
                                     contrast = contrast,
                                     geneList = geneList,
                                     geneListcolours = geneListcolours,
                                     customBg = customBg,
                                     comparisons = comparisons,
                                     plotOut = plotOut,
                                     plotType = plotType,
                                     pdfName = pdfName)
          resultsOut[[paste('codUp','comp',j,sep='_')]] <- selCodonOutUp
        }
        if(length(featSelDown) > 0){
          selCodonOutDown <- codonCalc(codonIn = codonIn,
                                   featsel = featSelDown,
                                   analysis = analysis,
                                   unit = unit,
                                   ads = ads,
                                   regulation = regulation,
                                   contrast = contrast,
                                   geneList = geneList,
                                   geneListcolours = geneListcolours,
                                   customBg = customBg,
                                   comparisons = comparisons,
                                   plotOut = plotOut,
                                   plotType = plotType,
                                   pdfName = pdfName)
          resultsOut[[paste('codDown','comp',j,sep='_')]] <- selCodonOutDown
        }
      }
    }
  }
  if(!is.null(addSign)){
    #
    signOutTmp <- signCalc(annot=annot,
                           addSign=addSign,
                           ads = ads,
                           geneList=geneList,
                           customBg=customBg)
    #
    resultsOut <- append(resultsOut, signOutTmp)
  }
  #
  for(cSel in unique(contrast)){
    contrastSel <- cSel
    
    genesSlopeFiltOut <- slopeFilt(ads=ads,
                                   regulationGen=regulationGen,
                                   contrastSel=contrastSel,
                                   minSlope=minSlope, 
                                   maxSlope=maxSlope)
                        
    miRNAOut <- miRNAanalysis(annot=annot,
                              miRNATargetScanFile = miRNATargetScanFile,
                              genesSlopeFiltOut = genesSlopeFiltOut,
                              ads=ads,
                              regulationGen=regulationGen,
                              contrastSel=contrastSel,
                              customBg=customBg,
                              rankIn=rankIn)
    #
    if(length(miRNAOut)>1){
      if(length(miRNAOut) > (miRNAmax)+1 ){
        miRNAOutFilt <- miRNAOut[2:(miRNAmax+1)]
      } else {
        miRNAOutFilt <- miRNAOut[2:length(miRNAOut)]
      }
      names(miRNAOutFilt) <- paste(names(miRNAOutFilt), paste('c',contrastSel,sep=''), sep='_')
      resultsOut <- append(resultsOut, miRNAOutFilt)
    }
  }
  
  ####Run feature integration
  utilsOut <- featureIntegration(features = resultsOut, 
                                 analysis_type = analysis_type,
                                 rmCat=rmCat,
                                 comparisons = comparisons,
                                 regOnly = regOnly, 
                                 allFeat = allFeat, 
                                 covarFilt = covarFilt,  
                                 NetModelSel= NetModelSel, 
                                 ads = ads, 
                                 regulationGen = regulationGen, 
                                 contrastSel = contrastSel, 
                                 geneList = geneList, 
                                 geneListcolours = geneListcolours, 
                                 customBg = customBg,
                                 effectMeasure = effectMeasure, 
                                 outDir = outDir,
                                 pdfName = pdfName)  

    return(utilsOut)
}




