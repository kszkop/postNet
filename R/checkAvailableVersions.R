
checkAvailableVersions <- function(species) {
  # Define the base directory for annotation files
  base_dir <- system.file("extdata/annotation/refseq", package = "postNet")
  
  # List existing species
  curr_tmp <- list.files(base_dir)
  
  # Check if the input species is a valid species name
  if (!species %in% curr_tmp) {
    stop("Invalid species name. The pre-prepared annotation file for that species does not exist. Please use a valid species name.")
  } else {
    # Create the directory path for the specified species
    species_dir <- file.path(base_dir, species)
    
    # List available versions for the specified species
    list_versions <- list.files(path = species_dir)
    
    if (length(list_versions) == 0) {
      cat("No annotation versions found for", species, "\n")
    } else {
      cat("Available versions for", species, ":\n")
      print(list_versions)
    }
  }
}