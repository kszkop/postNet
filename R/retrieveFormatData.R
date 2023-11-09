retrieveFormatData <- function(source,
                               version = NULL,
                               species = NULL,
                               customFile = NULL,
                               fastaFile = NULL,
                               posFile = NULL,
                               rna_gbff_file = NULL,
                               rna_fa_file = NULL,
                               genomic_gff_file = NULL) {
  
  # Validate the source input
  tryCatch({
    # Code that may throw an error
    checkSource(source)
  }, error = function(e) {
    stop("Source check failed: ", e$message)
  })
  
  # Validate specific parameters
  tryCatch({
    # Code that may throw an error
    checkInput(source, customFile, rna_gbff_file, rna_fa_file, genomic_gff_file, posFile)
  }, error = function(e) {
    stop("Parameters check failed: ", e$message)
  })

  # Check available species for the 'create' source
  if (source == "create") {
    if(!is_valid_species(species)){
      stop("Please specify a species, at the moment only 'human' or 'mouse' are available).") 
    }
    # Define URLs for downloading files
    url <- switch(species,
                  "human" = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/",
                  "mouse" = "https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/current/GCF_000001635.27-RS_2023_04/"
    )
    
    # Download files and unzip them
    download_and_unzip("customFasta.fa.gz", "GRCh38_latest_rna.fna.gz")
    download_and_unzip("customAnnot.gbff.gz", "GRCh38_latest_rna.gbff.gz")
    download_and_unzip("GeneRef.gff.gz", "GRCh38_latest_genomic.gff.gz")
    
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
    currTmp <- list.files(system.file("extdata/annotation/refseq", package = "anota2seqUtils"))
    
    if (!species %in% currTmp) {
      stop("This option is only available for species: human and mouse at the moment. Please use option createFromFile")
    }
    
    if (is.null(version)) {
      version <- checkAvailableVersions(species = species)
      versionInd <- sub("^[^.]*.", "", version)
      versionInd <- sort(versionInd, decreasing = TRUE)[1]
      version <- version[grep(versionInd, version)]
    }
    
    if (species == "human") {
      outDB <- read.delim(system.file(paste("extdata/annotation/refseq/human", version, sep = "/"), "humanDB.txt.gz", package = "anota2seqUtils"), stringsAsFactors = FALSE)
    }
    
    if (species == "mouse") {
      outDB <- read.delim(system.file(paste("extdata/annotation/refseq/mouse", version, sep = "/"), "mouseDB.txt.gz", package = "anota2seqUtils"), stringsAsFactors = FALSE)
    }
  } else if (source == "custom") {
    outDB <- read.delim(customFile, stringsAsFactors = FALSE)
    colnames(outDB) <- c('id', 'geneID', 'UTR5_seq', 'CDS_seq', 'UTR3_seq')
  } else if (source == "createFromFiles") {
    if(!is_valid_species(species)){
      stop("Please specify a species, at the moment only 'human' or 'mouse' are available).") 
    }
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
      url <- switch(species,
                    "human" = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/",
                    "mouse" = "https://ftp.ncbi.nlm.gov/refseq/M_musculus/annotation_releases/current/GCF_000001635.27-RS_2023_04/"
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
    stop("No correct option for the annotation file provided")
  }
  return(outDB)
}

