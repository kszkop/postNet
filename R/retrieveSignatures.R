#Prepare sequences and annotation
retrieveSignatures <- function(species=NULL #human and mouse only, it is required for source: create and load
){
  #
  #list existing species
  currTmp <- list.files(system.file("extdata/annotation/refseq",package = "anota2seqUtils"))
  #
  if(!species %in% currTmp){
    stop("This option is only  available for species: human and mouse at the moment")
  }
  if(species=='human'){
    load("../../jerryDexa/fromChristian/mRNA_signatures.Rdata")
    
    #load(system.file("extdata/annotation/refseq/human/", "mRNAsignatures_human.RData", package = "anota2seqUtils"))
  }
  if(species=='mouse'){
    #outDB <- read.delim(system.file("extdata/annotation/refseq/mouse", "mRNAsignatures_mouse.RData", package = "anota2seqUtils"), stringsAsFactors=FALSE)
  }
}