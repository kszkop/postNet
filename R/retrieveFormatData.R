retrieveFormatData <- function(source,
                               version = NULL,
                               species = NULL,
                               customFile = NULL,
                               fastaFile = NULL,
                               posFile = NULL,
                               rna_gbff_file = NULL,
                               rna_fa_file = NULL,
                               genomic_gff_file = NULL) {
  if (source == "create") {
    ####Check available species
    currTmp <- list.files(system.file("extdata/annotation/refseq/", package = "anota2seqUtils"))
    if (!species %in% currTmp) {
      stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
    }

    # Try to find solution so it looks automatically for the files. maybe sth with RCurl. Looks good just fomratting
    # if(species=='human'){
    # library("RCurl")
    # result <- getURL("https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/",verbose=TRUE,ftp.use.epsv=FALSE, dirlistonly = TRUE, crlf = TRUE)
    # result2 <- paste("https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/", strsplit(results, "\r*\n")[[1]], sep = "")
    # }
    # if(species=='mouse'){
    # Find current release
    # Path to ftp refSeq db
    
    ####Download files
    if (species == "human"
    ) {
      url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/"
      #
      download.file(paste(url, "GRCh38_latest_rna.fna.gz", sep = ""), destfile = "customFasta.fa.gz")
      download.file(paste(url, "GRCh38_latest_rna.gbff.gz", sep = ""), destfile = "customAnnot.gbff.gz")
      download.file(paste(url, "GRCh38_latest_genomic.gff.gz", sep = ""), destfile = "GeneRef.gff.gz")
    }
    if (species == "mouse") {
      url <- "https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/current/GCF_000001635.27-RS_2023_04/"
      #
      download.file(paste(url, "GCF_000001635.27_GRCm39_rna.fna.gz", sep = ""), destfile = "customFasta.fa.gz")
      download.file(paste(url, "GCF_000001635.27_GRCm39_rna.gbff.gz", sep = ""), destfile = "customAnnot.gbff.gz")
      download.file(paste(url, "GCF_000001635.27_GRCm39_genomic.gff.gz", sep = ""), destfile = "GeneRef.gff.gz")
    }
    R.utils::gunzip("customAnnot.gbff.gz")
    R.utils::gunzip("customFasta.fa.gz")
    R.utils::gunzip("GeneRef.gff.gz")

    ####Reformat
    seqs <- seqinr::read.fasta(file = "customFasta.fa", seqtype = "AA")
    seqs <- data.frame(id = sub("\\..*", "", names(seqs)), seq = t(as.data.frame(lapply(seqs, function(x) paste(x, collapse = "")))), row.names = NULL, stringsAsFactors = FALSE)
    
    if (species == "human") {
      command <- paste("perl", paste(perl.dir, "AnnotFromgbff_human.pl", sep = " "), intern = TRUE)
    }
    if (species == "mouse") {
      command <- paste("perl", paste(perl.dir, "AnnotFromgbff_mouse.pl", sep = " "), intern = TRUE)
    }
    system(command)
    
    annot <- read.delim("customAnnot.txt", stringsAsFactors = FALSE)
    colnames(annot) <- c("id", "UTR5_len", "CDS_stop", "Total_len")
    
    annotSeq <- merge(annot, seqs, by = "id")
    annotSeq <- extractRegSeq(annotSeq)

    gff <- gffRead("GeneRef.gff")
    bed <- extGff(gff)

    outDB <- merge(bed[, c(1, 7)], annotSeq[, c(1, 6, 7, 8)], by = "id")
    outDB <- outDB[grepl("NM_", outDB$id), ]

    write.table(outDB, file = "customDB.txt", col.names = T, row.names = F, sep = "\t", quote = F)
  } else if (source == "createFromSourceFiles") {
    # Unzip
    R.utils::gunzip(rna_gbff_file,remove=F)
    file.rename(gsub('.gz','',rna_gbff_file), 'customAnnot.gbff')
    R.utils::gunzip(rna_fa_file,remove=F)
    file.rename(gsub('.gz','',rna_fa_file), 'customFasta.fa')
    R.utils::gunzip(genomic_gff_file,remove=F)
    file.rename(gsub('.gz','',genomic_gff_file), 'GeneRef.gff')

    ####Reformat
    seqs <- seqinr::read.fasta(file = gsub('.gz','',rna_fa_file), seqtype = "AA")
    seqs <- data.frame(id = sub("\\..*", "", names(seqs)), seq = t(as.data.frame(lapply(seqs, function(x) paste(x, collapse = "")))), row.names = NULL, stringsAsFactors = FALSE)

    # if does not work yet. Change pl script so it takes any species and custom gbff file
    ##
    command <- paste("perl", paste(perl.dir, "AnnotFromgbff_human.pl", sep = " "), intern = TRUE)
    #
    system(command)
    #
    annot <- read.delim("customAnnot.txt", stringsAsFactors = FALSE)
    colnames(annot) <- c("id", "UTR5_len", "CDS_stop", "Total_len")
    #
    annotSeq <- merge(annot, seqs, by = "id")
    # Extract sequences:
    annotSeq <- extractRegSeq(annotSeq)

    # Extract info from gff
    gff <- gffRead(gsub('.gz','',genomic_gff_file))
    # Extract info from gff
    bed <- extGff(gff)

    outDB <- merge(bed[, c(1, 7)], annotSeq[, c(1, 6, 7, 8)], by = "id")
    outDB <- outDB[grepl("NM_", outDB$id), ]
    write.table(outDB, file = "customDB.txt", col.names = T, row.names = F, sep = "\t", quote = F)
  } else if (source == "load") {
    # list existing species
    currTmp <- list.files(system.file("extdata/annotation/refseq", package = "anota2seqUtils"))

    if (!species %in% currTmp) {
      stop("This option is only  available for species: human and mouse at the moment. Please use option createFromFile")
    }
    #
    if(is.null(version)){
      version <- checkAvailableVersions(species=species)
      #extract the latest
      versionInd <- sub("^[^.]*.","", version)
      versionInd <- sort(versionInd,decreasing = T)[1]
      version <- version[grep(versionInd, version)]
    }
    #
    if (species == "human") {
      outDB <- read.delim(system.file(paste("extdata/annotation/refseq/human", version, sep = "/"), "humanDB.txt.gz", package = "anota2seqUtils"), stringsAsFactors = FALSE)
    }
    if (species == "mouse") {
      outDB <- read.delim(system.file(paste("extdata/annotation/refseq/mouse", version, sep = "/"), "mouseDB.txt.gz", package = "anota2seqUtils"), stringsAsFactors = FALSE) # }
    }
  } else if (source == "custom") {
    outDB <- read.delim(customFile, stringsAsFactors = FALSE)
    colnames(outDB) <- c('id', 'geneID', 'UTR5_seq', 'CDS_seq', 'UTR3_seq')
  } else if (source == "createFromFiles") {
    # Load files
    posTmp <- read.delim(posFile, stringsAsFactors = FALSE)
    colnames(posTmp) <- c("id", "UTR5_len", "CDS_stop", "Total_len")
    #
    seqs <- seqinr::read.fasta(fastaFile, seqtype = "AA")
    seqs <- data.frame(id = sub("\\..*", "", names(seqs)), seq = t(as.data.frame(lapply(seqs, function(x) paste(x, collapse = "")))), row.names = NULL, stringsAsFactors = FALSE)
    #
    annotSeq <- merge(posTmp, seqs, by = "id")
    # Extract sequences:
    annotSeq <- extractRegSeq(annotSeq)
    #
    # Path to ftp refSeq db
    if (is.null(genomic_gff_file)) {
      #
      if (species == "human") {
        url <- "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/"
        download.file(paste(url, "GRCh38_latest_genomic.gff.gz", sep = ""), destfile = "GeneRef.gff.gz")
      }
      if (species == "mouse") {
        url <- "https://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/annotation_releases/current/GCF_000001635.27-RS_2023_04/"
        download.file(paste(url, "GCF_000001635.27_GRCm39_genomic.gff.gz", sep = ""), destfile = "GeneRef.gff.gz")
      }
      R.utils::gunzip("GeneRef.gff.gz")
      #
      gff <- gffRead("GeneRef.gff")
    } else {
      gff <- gffRead(genomic_gff_file)
    }
    # Extract info from gff
    bed <- extGff(gff)
    #
    outDB <- merge(bed[, c(1, 7)], annotSeq[, c(1, 6, 7, 8)], by = "id")
    outDB <- outDB[grepl("NM_", outDB$id), ]
    write.table(outDB, file = "customDB.txt", col.names = T, row.names = F, sep = "\t", quote = F)
  } else {
    stop("No correct option for annotation file provided")
  }
  return(outDB)
}
