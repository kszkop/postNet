postNetStart <- function(ads = NULL,
                                regulation = NULL,
                                contrast = NULL,
                                regulationGen = NULL, 
                                contrastSel = NULL,
                                geneList = NULL,
                                geneListcolours = NULL,
                                customBg = NULL,
                                effectMeasure = NULL,
                                selection = "random",
                                source,
                                version = NULL,
                                species = NULL,
                                customFile = NULL,
                                fastaFile = NULL,
                                posFile = NULL,
                                rna_gbff_file = NULL,
                                rna_fa_file = NULL,
                                genomic_gff_file = NULL,
                                adjObj = NULL,
                                region_adj = NULL,
                                excl = FALSE,
                                keepAll = FALSE) {

  #
  if(!is.null(ads) && !is.null(geneList)){
    stop("Please provide an anota2seq object or a genelist, not both.")
  }
  if(!is.null(ads)){
    check_ads(ads)
    
    if (!is.null(regulation) && !is.character(regulation) && !regulation %in% c("translationUp","translationDown","translatedmRNAUp","translatedmRNADown","bufferingmRNAUp","bufferingmRNADown","mRNAAbundanceUp","mRNAAbundanceDown","totalmRNAUp","totalmRNADown")) {
      stop("'regulation' should be a character vector chosen from translationUp,translationDown,translatedmRNAUp,translatedmRNADown,bufferingmRNAUp,bufferingmRNADown,mRNAAbundanceUp,mRNAAbundanceDown,totalmRNAUp,totalmRNADown")
    }
    if (!is.null(regulation)){
      if(!is.null(contrast) && !is.numeric(contrast) && !length(contrast) == length(regulation) && !contrast %in% seq(1,ncol(ads@contrasts),1)){
        stop("The 'contrast' provided should be a numeric vector corresponding to the numbers of the anota2seq contrasts for each regulatory mode selected for comparison. Please see the anota2seq vignette for additional details on contrasts.")
      }
    }
  }
  if(is.null(ads)){
    if(is.null(geneList)){
      stop('Please provide either an anota2seq object or a gene list.')
    } else {
      check_geneList(geneList)
      
      if (!is.null(geneListcolours) && !is.character(geneListcolours) && !length(geneListcolours)== length(geneList)) {
        stop("'geneListcolours' should be a character vector of the same length as geneList. These colours will be used for plotting.")
      }
    }
  }
  if(!is.null(customBg)){
    if(!is.character(customBg)){
      stop("The 'customBg' provided must be a character vector.")
    }
    if(!length(setdiff(unlist(geneList), customBg))==0){
      stop("There are entries in geneList that are not in 'customBg'. Please ensure that all genes are included in the background.")
    }
  }
  check_source(source)
  check_selection(selection)
  check_input(source, customFile, rna_gbff_file, rna_fa_file, genomic_gff_file, posFile)
  
  if(!is.null(adjObj)){
    check_adjObj(adjObj)
    valid_regions <- c('UTR5', 'UTR3')
    if (!all(region_adj %in% valid_regions)) {
      stop("'region_adj' must be provided and can only inlcude: 'UTR5', and/or 'UTR3'. It should also match the named entries in the list adjObj.")
    }
    if(!is_logical(excl)){
      stop("'excl' must be logical: TRUE of FALSE.")
    }
    if(!is_logical(keepAll)){
      stop("'keepAll' must be logical: TRUE of FALSE.")
    }
  }
  # Check available species for the 'create' source
  if (source == "create") {
    if(!is_valid_species(species)){
      stop("Please specify a species, at the moment only 'human' or 'mouse' are available).") 
    }
    ####Download files
    if (species == "human"
    ) {
      url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/current/"
      
      responseTmp <- httr::GET(url)
      if (httr::status_code(responseTmp) == 200) {
        
        pageTmp <- httr::content(responseTmp, as = "text") %>% rvest::read_html()
        linksTmp <- pageTmp %>% rvest::html_nodes("a") %>% rvest::html_attr("href")
        #
        version <- linksTmp[grepl("GCF_000001405", linksTmp)]
      } else {
        cat("Failed to fetch the URL. Status code:", httr::status_code(responseTmp), "\n")
      }
      #
      
      download.file(paste(url, version, "_rna.fna.gz", sep = ""), destfile = "customFasta.fa.gz")
      download.file(paste(url, version, "_rna.gbff.gz", sep = ""), destfile = "customAnnot.gbff.gz")
      download.file(paste(url, version, "_genomic.gff.gz", sep = ""), destfile = "GeneRef.gff.gz")
    }
    if (species == "mouse") {
      url <- "https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/current/"
      
      responseTmp <- httr::GET(url)
      if (httr::status_code(responseTmp) == 200) {

        pageTmp <- httr::content(responseTmp, as = "text") %>% rvest::read_html()
        linksTmp <- pageTmp %>% rvest::html_nodes("a") %>% rvest::html_attr("href")
        #
        version <- linksTmp[grepl("GCF_000001635", linksTmp)]
      } else {
        cat("Failed to fetch the URL. Status code:", httr::status_code(responseTmp), "\n")
      }
      #
      download.file(paste(url, version,"_rna.fna.gz", sep = ""), destfile = "customFasta.fa.gz")
      download.file(paste(url, version,"_rna.gbff.gz", sep = ""), destfile = "customAnnot.gbff.gz")
      download.file(paste(url, version,"_genomic.gff.gz", sep = ""), destfile = "GeneRef.gff.gz")
    }
    R.utils::gunzip("customAnnot.gbff.gz")
    R.utils::gunzip("customFasta.fa.gz")
    R.utils::gunzip("GeneRef.gff.gz")

    # Reformat sequence data
    seqs <- seqinr::read.fasta(file = "customFasta.fa", seqtype = "AA")
    seqs <- data.frame(id = sub("\\..*", "", names(seqs)),
                       seq = t(as.data.frame(lapply(seqs, function(x) paste(x, collapse = "")))),
                       row.names = NULL,
                       stringsAsFactors = FALSE)
    
    # Determine the appropriate Perl script based on species
    perl_script <- switch(species,
                          "human" = "AnnotFromgbff_human.pl",
                          "mouse" = "AnnotFromgbff_mouse.pl"
    )
    # Run the Perl script
    command <- paste("perl", file.path(perl.dir, perl_script), sep = " ")
    system(command)
    
    # Read and merge annotation data
    annot <- read.delim("customAnnot.txt", stringsAsFactors = FALSE)
    colnames(annot) <- c("id", "UTR5_len", "CDS_stop", "Total_len")
    
    annotSeq <- merge(annot, seqs, by = "id")
    annotSeq <- extractRegSeq(annotSeq)
    
    gff <- gffRead("GeneRef.gff")
    bed <- extGff(gff)
    
    outDB <- merge(bed[, c(1, 7)], annotSeq[, c(1, 6, 7, 8)], by = "id")
    outDB <- outDB[grepl("NM_", outDB$id), ]
    
    write.table(outDB, file = "customDB.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  } else if (source == "createFromSourceFiles") {
    if(!is_valid_species(species)){
      stop("Please specify a species, at the moment only 'human' or 'mouse' are available).") 
    }
    # Unzip the source files
    source_files <- c(rna_gbff_file, rna_fa_file, genomic_gff_file)
    source_files <- gsub('.gz', '', source_files)
    filenames <- c("customAnnot.gbff", "customFasta.fa", "GeneRef.gff")
    for (i in 1:length(source_files)) {
      R.utils::gunzip(source_files[i], remove = FALSE)
      file.rename(source_files[i], filenames[i])
    }
    
    # Reformat sequence data
    seqs <- seqinr::read.fasta(file = "customFasta.fa", seqtype = "AA")
    seqs <- data.frame(id = sub("\\..*", "", names(seqs)),
                       seq = t(as.data.frame(lapply(seqs, function(x) paste(x, collapse = "")))),
                       row.names = NULL,
                       stringsAsFactors = FALSE)
    
    # Run the Perl script (adjust as needed)
    # Determine the appropriate Perl script based on species
    perl_script <- switch(species,
                          "human" = "AnnotFromgbff_human.pl",
                          "mouse" = "AnnotFromgbff_mouse.pl"
    )
    # Run the Perl script
    command <- paste("perl", file.path(perl.dir, perl_script), sep = " ")
    system(command)
    
    # Read and merge annotation data
    annot <- read.delim("customAnnot.txt", stringsAsFactors = FALSE)
    colnames(annot) <- c("id", "UTR5_len", "CDS_stop", "Total_len")
    
    annotSeq <- merge(annot, seqs, by = "id")
    annotSeq <- extractRegSeq(annotSeq)
    
    # Extract info from the GFF file
    gff <- gffRead("GeneRef.gff")
    bed <- extGff(gff)
    
    outDB <- merge(bed[, c(1, 7)], annotSeq[, c(1, 6, 7, 8)], by = "id")
    outDB <- outDB[grepl("NM_", outDB$id), ]
    
    write.table(outDB, file = "customDB.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  } else if (source == "load") {
    if(!is_valid_species(species)){
      stop("Please specify a species, at the moment only 'human' or 'mouse' are available).") 
    }
    # List existing species
    currTmp <- list.files(system.file("extdata/annotation/refseq", package = "postNet"))
    
    if (!species %in% currTmp) {
      stop("This option is currently only available for species 'human' and 'mouse'. Please use the options 'custom' and 'customFile' to provide annotations for other species.")
    }
    
    if (is.null(version)) {
      version <- checkAvailableVersions(species = species)
      versionInd <- sub("^[^.]*.", "", version)
      versionInd <- sort(versionInd, decreasing = TRUE)[1]
      version <- version[grep(versionInd, version)]
    }
    if (species == "human") {
      outDB <- read.delim(system.file(paste("extdata/annotation/refseq/human", version, sep = "/"), "humanDB.txt.gz", package = "postNet"), stringsAsFactors = FALSE)
    }
    if (species == "mouse") {
      outDB <- read.delim(system.file(paste("extdata/annotation/refseq/mouse", version, sep = "/"), "mouseDB.txt.gz", package = "postNet"), stringsAsFactors = FALSE)
    }
  } else if (source == "custom") {
    outDB <- read.delim(customFile, stringsAsFactors = FALSE)
    colnames(outDB) <- c('id', 'geneID', 'UTR5_seq', 'CDS_seq', 'UTR3_seq')
  } else if (source == "createFromFasta") {
    
    posTmp <- read.delim(posFile, stringsAsFactors = FALSE)
    colnames(posTmp) <- c("id", "UTR5_len", "CDS_stop", "Total_len")
    
    seqs <- seqinr::read.fasta(fastaFile, seqtype = "AA")
    seqs <- data.frame(id = sub("\\..*", "", names(seqs)),
                       seq = t(as.data.frame(lapply(seqs, function(x) paste(x, collapse = "")))),
                       row.names = NULL,
                       stringsAsFactors = FALSE)
    
    annotSeq <- merge(posTmp, seqs, by = "id")
    annotSeq <- extractRegSeq(annotSeq)
    
    # Path to ftp refSeq db
    if (is.null(genomic_gff_file)) {
      if(!is_valid_species(species)){
        stop("Please specify a species. Currently, 'human' or 'mouse' are available. For use with other species annotations, please use the 'custom' option and provide an annotation file using the 'customeFile' argument.") 
      }
      #
      url <- switch(species,
                    "human" = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/",
                    "mouse" = "https://ftp.ncbi.nlm.gov/refseq/M_musculus/annotation_releases/current/GCF_000001635.27-RS_2024_02/"
      )
      
      download.file(paste(url, switch(species, "human" = "GRCh38_latest_genomic.gff.gz", "mouse" = "GCF_000001635.27_GRCm39_genomic.gff.gz"), sep = ""), destfile = "GeneRef.gff.gz")
      R.utils::gunzip("GeneRef.gff.gz")
      
      gff <- gffRead("GeneRef.gff")
    } else {
      gff <- gffRead(genomic_gff_file)
    }
    
    bed <- extGff(gff)
    
    outDB <- merge(bed[, c(1, 7)], annotSeq[, c(1, 6, 7, 8)], by = "id")
    outDB <- outDB[grepl("NM_", outDB$id), ]
    
    write.table(outDB, file = "customDB.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  } else {
    stop("Please select one of: 'custom', 'create', 'createFromFasta', and 'createFromSourceFiles' for the 'source' argument. See the help page and vignette for additional details on sequence annotations.")
  }
  
  ####
  annotBg <- gSel(annot = outDB, ads = ads, customBg = customBg, geneList = geneList)
  #
  if(!is.null(adjObj)){
    annotBg <- adjustSeq(annot=annotBg, region_adj = region_adj, adjObj = adjObj, keepAll =  keepAll, excl = excl)
  }
  
  annot <- new("postNetAnnot",
              UTR5 = NULL,
              CDS = NULL,
              UTR3 = NULL,
              CCDS = NULL)
  
  for(reg in c('UTR5','CDS','UTR3')){
    annotTmp <- regSel(annot = annotBg, region = reg)
    annotBgSelTmp <- isoSel(annot = annotTmp, method = selection)
    
    RegionTmp <- new("postNetRegion",
                     id = annotBgSelTmp$id,
                     geneID = annotBgSelTmp$geneID,
                     sequences = annotBgSelTmp$seqTmp)
    
    if(reg == 'UTR5'){
      annot@UTR5 <- RegionTmp
    } else if (reg == 'CDS'){
      annot@CDS <- RegionTmp
    } else if (reg == 'UTR3'){
      annot@UTR3 <- RegionTmp
    }
  }
  
  ##
  genesIn <- resSel(ads = ads, regulation = regulation, contrast = contrast, geneList = geneList)
  
  #add here to check numbers of genes
  #if(length(resOut)==0){
  #  stop('There are no regulated genes. Check the input or run without indicating regulation and comparisons')
  #}
  
  coloursIn <- coloursSel(ads = ads, genesIn = genesIn, geneList = geneList, geneListcolours = geneListcolours)
  effIn <- effectSel(ads = ads, regulationGen = regulationGen, contrastSel = contrastSel, effectMeasure = effectMeasure)
  bgIn <- getBg(ads = ads, customBg = customBg, geneList = geneList)
  #
  dataIn <- new("postNetDataIn",
                background = bgIn,
                geneList = genesIn,
                effect = effIn,
                colours = coloursIn)
  
  analysis <- new("postNetAnalysis",
                  featureIntegration = NULL,
                  motifs= NULL,
                  codons = NULL,
                  GO = NULL,
                  GSEA = NULL,
                  GAGE = NULL,
                  miRNA = NULL)
  
  # initialize the postNet
  postNetData <- new("postNetData",
                             version = version,
                             species = species,
                             selection = selection,
                             annot = annot,
                             dataIn = dataIn,
                             features = NULL,
                             analysis = analysis)
  
  
  return(postNetData)
}

