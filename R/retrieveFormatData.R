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
  
  # Validate species input
  tryCatch({
    # Code that may throw an error
    checkSpecies(source, species)
  }, error = function(e) {
    stop("Species check failed: ", e$message)
  })
  # Validate specific parameters
  tryCatch({
    # Code that may throw an error
    checkParameters(source, customFile, rna_gbff_file, rna_fa_file, genomic_gff_file, posFile)
  }, error = function(e) {
    stop("Parameters check failed: ", e$message)
  })

  # Check available species for the 'create' source
  if (source == "create") {
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
    command <- paste("perl", file.path(perl.dir, "AnnotFromgbff_human.pl"), sep = " ")
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

# Function to validate the source input
checkSource <- function(source) {
  valid_sources <- c("create", "createFromSourceFiles", "load", "custom", "createFromFiles")
  if (!(source %in% valid_sources)) {
    stop("Invalid source. Please provide a valid source option.")
  }
}

# Function to validate the species input
checkSpecies <- function(source, species) {
  if (source %in c("create", "load") && is.null(species)) {
    stop("Please specify a species (e.g., 'human' or 'mouse').")
  }
  if (source %in c("create", "load") && !(species %in% c("human", "mouse"))) {
    stop("This option is only available for species: human and mouse at the moment. Please use option createFromFile.")
  }
}

# Function to validate specific input parameters
checkParameters <- function(source, customFile, rna_gbff_file, rna_fa_file, genomic_gff_file, posFile) {
  if (source == "createFromSourceFiles") {
    if (is.null(rna_gbff_file)) {
      stop("Please provide an rna_gbff_file.")
    }
    if (is.null(rna_fa_file)) {
      stop("Please provide an rna_fa_file.")
    }
    if (is.null(genomic_gff_file)) {
      stop("Please provide a genomic_gff_file.")
    }
  } else if (source == "custom") {
    if (is.null(customFile)) {
      stop("Please provide a customFile.")
    }
  } else if (source == "createFromFiles") {
    if (is.null(posFile)) {
      stop("Please provide a posFile in the format: id, UTR5_len, CDS_stop, Total_len.")
    }
  }
}

# Function to download and unzip files
downloadAndUnzip <- function(source, species) {
  url <- generateDownloadURL(species)
  download_and_unzip("customFasta.fa.gz", paste0(url, "_rna.fna.gz"))
  download_and_unzip("customAnnot.gbff.gz", paste0(url, "_rna.gbff.gz"))
  download_and_unzip("GeneRef.gff.gz", paste0(url, "_genomic.gff.gz"))
}
