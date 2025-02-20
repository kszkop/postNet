\name{postNetStart}
\alias{postNetStart}
\alias{postNetData-class}
\title{Retrieve and format reference sequences and inputs for postNet analysis
}
\description{
  The \code{postNetStart} function is the first step in a \pkg{postNet} analysis workflow. It retrieves and compiles    all necessary inputs and annotations into a \code{postNetData} S4 object, including:
    \itemize{
      \item Reference sequence annotations (5'UTR, CDS, 3'UTR), including in-built, custom, or retrieved directly 
      from the NCBI RefSeq database. Optionally, UTR sequences can also be adjusted if more precise sequences are 
      available (e.g. from CAGE, QuantSeq, etc.)
      \item Gene sets of interest (either from an \pkg{anota2seq} analysis by providing the \code{Anota2seqDataSet}
      object, or using custom gene lists of interest).
      \item A numeric regulation effect measurement (e.g., fold changes in gene expression or translation 
      efficiency, etc.) either from an \pkg{anota2seq} analysis, or from a custom source. These values will be used
      in the \link[postNet]{featureIntegration} analysis, where mRNA features, gene signatures, and/or other variables
      can be used to explain changes in the regulation effect measurement using step-wise regression and/or random
      forest approaches.
    }
  The resulting \code{postNetData} object is then passed to downstream functions in the \pkg{postNet}
  workflow, and stores analysis results.
}
\usage{
  postNetStart(
    ## ----------------------------
    ## 1. Input Data (Anota2seq or Custom)
    ## ----------------------------
    ads = NULL,
    regulation = NULL,
    contrast = NULL,
    regulationGen = NULL,
    contrastSel = NULL,
    geneList = NULL,
    geneListcolours = NULL,
    customBg = NULL,
    effectMeasure = NULL,
    
    ## ----------------------------
    ## 2. Reference Sequence Annotations
    ## ----------------------------
    source,
    species = NULL,
    version = NULL,
    customFile = NULL,
    fastaFile = NULL,
    posFile = NULL,
    rna_gbff_file = NULL,
    rna_fa_file = NULL,
    genomic_gff_file = NULL,
    selection = "random",
    setSeed = NULL,
    
    ## ----------------------------
    ## 3. Adjusting Sequences (Custom UTRs)
    ## ----------------------------
    adjObj = NULL,
    region_adj = NULL,
    excl = FALSE,
    keepAll = FALSE
  )
}

\arguments{
  \item{ads}{An optional S4 object of class \code{Anota2seqDataSet} (resulting from the \code{\link[anota2seq]{anota2seqRun}} function in the \pkg{anota2seq} package). Lists of regulated genes (regulatory modes) and regulatory effect measures will be retrieved from this object for the analysis. \strong{Do not} provide both \code{ads} and \code{geneList} at the same time.
}
  \item{regulation}{If using \code{ads}, provide a character vector of regulatory modes to retrieve from the \code{Anota2seqDataSet} object. Valid options include any individual selection, or combination of: \code{"translationUp"}, \code{"translationDown"}, \code{"translatedmRNAUp"}, \code{"translatedmRNADown"}, \code{"bufferingmRNAUp"}, \code{"bufferingmRNADown"}, \code{"mRNAAbundanceUp"}, \code{"mRNAAbundanceDown"}, \code{"totalmRNAUp"}, \code{"totalmRNADown"}. Please see the \pkg{anota2seq} vignette for additional details on regulatory modes.}
  \item{contrast}{If using \code{ads}, provide a numeric vector specifying which \pkg{anota2seq} contrast each element in \code{regulation} should be selected from. For example, \code{c(1,1,2)} selects the first two regulatory modes provided in \code{regulation} from contrast 1, and the third mode from contrast 2. This vector must always be the same length as the input for \code{regulation}, even if the \code{Anota2seqDataSet} has only one contrast. For more details on contrasts, please see the \pkg{anota2seq} vignette.
}
  \item{regulationGen}{If using \code{ads}, select the "general" regulation effect measurement that will be taken from the \code{ads} object. These values will be used in the downstream feature integration where the mRNA features (and/or other input signatures) quantified in the specified gene sets (\code{regulation}) will be used to explain changes in the regulatory effect measurement (i.e. translation efficiency, etc.) with the \code{\link[postNet]{featureIntegration}} function. This parameter can be \strong{one} of: \code{"translation"}, \code{"buffering"}, \code{"mRNAAbundance"}, \code{"translatedmRNA"}, or \code{"totalmRNA"}.
}
  \item{contrastSel}{If using \code{ads}, select the \pkg{anota2seq} contrast where the "general" regulation effect measurement corresponding to \code{regulationGen} should be extracted from. This must be a single numeric value. 
}
  \item{geneList}{If \pkg{anota2seq} was not used, instead of providing \code{ads}, gene sets of interest can be provided as a named list of one or more character vectors. Each list element should be a set of genes (For example, \code{list(set1=c("GeneA","GeneB"), set2=c("GeneC","GeneD"))}). Ensure that the gene IDs provided match those in the reference sequence annotations.  
}
  \item{geneListcolours}{If using \code{geneList}, provide a character vector of colors to be used for plotting. This must have the same length as \code{geneList}. If using \code{ads}, the colours for each regulatory mode from \pkg{anota2seq} will be used automatically. 
}
  \item{customBg}{If using \code{geneList}, optionally provide a character vector of gene/transcript IDs corresponding to the custom background. This background will be used for statistical comparisons against the gene sets of interest provided in \code{geneList}. Note that the background should correspond to the set of genes or transcripts that was quantified (i.e., all genes used as input in the analysis that generated the regulation effect measurement [\code{effectMeasure}]). For example, all genes in a given data set passing expression/reproducibility thresholds. Note that all genes in \code{geneList} must be found in \code{customBg}. If \code{customBg} is \code{NULL}, the background will be set as all genes IDs present in the provided \code{geneList}. Although optional, it is strongly recommended to carefully consider that the appropriate background is used. If using \code{ads}, the appropriate background is automatically retrieved from the \code{Anota2seqDataSet} object, and \code{customBg} should be \code{NULL}. 
}
  \item{effectMeasure}{If using \code{geneList}, provide a numeric vector corresponding to a custom regulation effect measurement. For example, fold changes in mRNA expression or translation efficiency (although any continuous numeric variable could theoretically be used). The elements of the vector should be named with the gene/transcript IDs in a format corresponding to those provided in the \code{geneList}. This parameter should not be provided when using \code{ads}, as the regulation effect measurement will be retrieved from the \code{Anota2seqDataSet} using the \code{regulationGen} parameter. 
}
  \item{selection}{Specifies which mRNA isoform will be selected from the reference sequence annotation for use in analyses in cases when multiple isoforms are present for a given gene. This parameter can be one of:
\describe{
     \item{\code{"longest"}}{Selects the longest sequence isoform.}
     \item{\code{"shortest"}}{Selects the shortest sequence isoform.}
     \item{\code{"random"}}{(Default) Selects an isoform at random. Note that selecting the longest or
     shortest isoforms for all genes may skew results for some analyses. For this reason, the default
     selection is \code{"random"}. However, this may lead to slight differences in results if the \code{postNetStart} 
     function is run multiple times. For cases where multiple analyses may be run using the same data, the \code{setSeed}
     parameter can be used to ensure that random isoform selection is stable each time \code{postNetStart} is run.}
     }
}
  \item{setSeed}{If \code{selection} is \code{"random"}, provide a single integer value to be used with \code{set.seed} 
  to ensure that the random isoform selection is reproducible between different runs of the \code{postNetStart} function. 
  Any integer can be provided (for example 123), however, the same value must be provided between different analyses for 
  isoform selection to be reproducible. 
}
  \item{source}{Select the source of the reference sequence annotations. Note that the the gene/transcript annotations used in the \code{ads} or \code{geneList} should be the same as those used in the reference sequence annotations. This parameter can be set to:
   \describe{
     \item{\code{"load"}}{(Default) Load a pre-made RefSeq reference sequence annotation of the selected 
         \code{species} and \code{version} included in \pkg{postNet}. Currently, this option is supported with
         the \code{species} \code{"human"} and \code{"mouse"}. Available RefSeq versions can be checked using the
         \code{\link[postNet]{checkAvailableVersions}} function.}
     \item{\code{"create"}}{Automatically download the most recent version of the RefSeq annotation files 
         from the NCBI database for the selected \code{species} and build a new reference sequence annotation locally.
         Currently, this option is supported with the \code{species} \code{"human"} and \code{"mouse"}.}
     \item{\code{"createFromSourceFiles"}}{Use local RefSeq reference files already downloaded from the NCBI RefSeq
         database to create annotation sequences. Three files are required, the RNA GBFF ("rna.gbff"), RNA FASTA ("rna
         .fna"), and genomic GFF ("genomic.gff"). These files must be provided using the \code{rna_gbff_file},
         \code{rna_fa_file}, and \code{genomic_gff_file} parameters. This option can be used to assemble reference
         sequence annotations for RefSeq versions not currently supported with \code{"load"}, and unlike the
         \code{"create"} option does not require an internet connection if the reference files have already been
         downloaded. Currently, this option is supported with the \code{species} \code{"human"} and \code{"mouse"}.}
     \item{\code{"createFromFasta"}}{Alternatively, reference sequence annotations can be assembled from a FASTA file 
         (provided with the \code{fastaFile} parameter) with reference sequences, and a file specifying the positions of 
         the sequence regions (provided with the \code{posFile} parameter). Optionally provide a genomic GFF file (with
         the \code{genomic_gff_file} parameter). If no GFF is provided the latest version for \code{species} will be 
         automatically downloaded from the NCBI RefSeq database as for \code{"create"}. Currently, this option is
         supported with the \code{species} \code{"human"} and \code{"mouse"}.}
     \item{\code{"custom"}}{Use a pre-prepared reference sequence annotation file provided with the \code{customFile}
         parameter. This option can be used when working with sequences from \code{species} not currently supported
         by the options above, or with annotations other than those in the NCBI RefSeq database.}
   }
}
  \item{version}{If \code{source} is \code{"load"}, specify the RefSeq release version that should be used for the reference sequence annotations. Available versions can be checked using the \code{\link[postNet]{checkAvailableVersions}}
function. If \code{NULL}, the most recent version is used.
}
  \item{species}{Select the species of interest. Currently, \code{"human"} and \code{"mouse"} are supported. This parameter is required unless \code{source} is \code{"custom"}.
}
  \item{customFile}{If \code{source} is \code{"custom"}, provide the name and file path to the custom reference sequence
annotation file. The \code{customFile} must be TAB-delimited, with the headings: \code{transcriptID, geneID, UTR5_seq, CDS_seq, UTR3_seq}. 
}
  \item{fastaFile}{If \code{source} is \code{"createFromFasta"}, provide the name and file path to a FASTA file to be used to assemble reference sequence annotations.
}
  \item{posFile}{If \code{source} is \code{"createFromFasta"}, provide the name and file path to a file specifying positions of transcript sequence regions in the \code{fastaFile}. The \code{posFile} must be TAB-delimited, with the headings: \code{transcriptID, utr5_stop, cds_stop, total_length}. Please see the vignette for additional details.
}
  \item{rna_gbff_file}{If \code{source} is \code{"createFromSourceFiles"}, provide the name and file path to a RNA GBFF reference file for either human or mouse. This can can be downloaded from the NCBI RefSeq database. Note that the file should be compressed with gzip (.gz extension), which is the default format for NCBI RefSeq downloads.
}
  \item{rna_fa_file}{If \code{source} is \code{"createFromSourceFiles"}, provide the name and file path to a RNA FASTA reference file for either human or mouse. This can can be downloaded from the NCBI RefSeq database. Note that the file should be compressed with gzip (.gz extension), which is the default format for NCBI RefSeq downloads.
}
  \item{genomic_gff_file}{If \code{source} is \code{"createFromSourceFiles"}, provide the name and file path to a genomic GFF reference file for either human or mouse. This can can be downloaded from the NCBI RefSeq database. Note that the file should be compressed with gzip (.gz extension). If \code{source} is \code{"createFromFasta"}, an unzipped genomic GFF can be specified manually using this parameter. Alternatively, leaving \code{NULL} will automatically download the most recent version for \code{species} from the NCBI RefSeq database.
}
  \item{adjObj}{A named list of custom UTR sequences that will replace the existing UTR sequences in the annotation file. List elements must be named \code{"UTR5"} and/or \code{"UTR3"}, and elements should be a named character vector of UTR sequences where names are transcript IDs (for example, \code{NM_...}), and values are the new UTR sequences. This parameter is useful when more precise measurements of UTR sequences are available (for example, from CAGE or QuantSeq data).
}
  \item{region_adj}{Character vector specifying which UTR regions in your annotation should be adjusted by the provided
\code{adjObj}. This parameter can be used to adjust both UTR regions (e.g., \code{c("UTR5","UTR3")}), or one at a time.
}
  \item{excl}{When providing custom UTR sequences with \code{adjObj}, sequences may not be available for all genes/isoforms included in the existing annotation file. This parameter determines whether to keep all genes/isoforms, or exclude those that are not included in \code{adjObj}. If \code{TRUE}, any transcript IDs present in the existing annotation file that are not present in \code{adjObj} will be removed. If \code{FALSE}, then the UTR sequences in the annotation file will be replaced with the custom sequences in \code{adjObj} when they are available, but for transcript IDs where a custom sequence is not available, the original sequence present in the annotation file will be retained and used in the analysis. The default is \code{FALSE}.
}
  \item{keepAll}{When using adjusting UTR sequences with \code{adjObj} and \code{excl} is \code{FALSE}, the \code{keepAll} parameter determines how the custom UTR sequences are integrated into the existing annotation file for genes with multiple isoforms. If \code{TRUE}, the UTR sequences in the annotation file will be replaced with the custom sequences corresponding to the same transcript IDs. However, if there are multiple transcript IDs (isoforms) for the same gene ID, all isoforms will be retained, including transcript IDs where custom sequences were not available in \code{adjObj}. Note that unless the custom sequences are always the longest or shortest isoforms that can be selected using the \code{selection} parameter, when \code{keepAll} is \code{TRUE}, the custom sequences are not guaranteed to be the ones considered if a gene-level analysis is performed. If \code{FALSE}, transcript IDs (isoforms) corresponding to the same gene that are not found in the custom sequences provided in \code{adjObj} will be removed, such that only the custom UTR sequences will be available for genes with multiple isoforms. The default is \code{FALSE}.
}
}
\details{
  By default, \code{postNetStart} will load or create reference sequence 
  annotations from RefSeq for \code{species}, assemble or load gene lists of 
  interest, define a suitable background gene set, and store a numeric regulation 
  effect measurement (e.g., fold change) for downstream modeling. 
  
  \strong{Reproducibility Note}: If \code{selection = "random"}, repeated calls 
  without \code{setSeed} can produce different isoform selections. To ensure 
  reproducible results between different analyses of the same dataset, set \code{setSeed} 
  to a fixed integer.
  
  When adjusting annotations via \code{adjObj}, ensure that transcript IDs 
  match those in your reference annotation. You can choose to exclude or keep 
  genes that lack updated sequences via \code{excl} and \code{keepAll}.
  
  The outputs of each step of postNet analysis will also be stored in the \code{postNetData}
  object, and can be retrieved through the use of a series of helper functions (see Note and examples below).
  
  See the \pkg{postNet} vignette for a more extensive tutorial and examples.
}

\value{An S4 object of class \code{postNetData} containing:
    \itemize{
      \item \strong{Reference Sequence Annotations:} 5'UTR, CDS, and 3'UTR sequences for each 
      gene and/or transcript isoform, with the possiblity to customize UTRs via \code{adjObj}.
      \item \strong{Gene Lists:} One or more sets of gene IDs of interest, as well as a set of
      gene IDs to be used as background for statistical comparisons (from \code{ads}, or 
      \code{geneList} and \code{customBg}).
      \item \strong{Regulation Effect Measurement:} A numeric vector of an effect of interest,
      (e.g., translation efficiency fold changes) that will be modelled using mRNA features and/or
      other signatures of interest.
      \item \strong{Metadata:} The chosen \code{species}, sequence annotation \code{version}, and
      isoform \code{selection} method.
      \item \strong{Analysis Results:} The results of subsequent steps of the
      \pkg{postNet} workflow, such as \code{\link[postNet]{featureIntegration}}, are stored in the
      \code{postNetData} object.
    }
}

\note{The helper functions available to retrieve inputs and results from the \code{postNetData} object are:
  \describe{
     \item{\link[postNet]{ptn_backround}}{}
     \item{\link[postNet]{ptn_codonAnalysis}}{}
     \item{\link[postNet]{ptn_codonSelection}}{}
     \item{\link[postNet]{ptn_colours}}{}
     \item{\link[postNet]{ptn_dataIn}}{}
     \item{\link[postNet]{ptn_effect}}{}
     \item{\link[postNet]{ptn_features}}{}
     \item{\link[postNet]{ptn_GAGE}}{}
     \item{\link[postNet]{ptn_geneID}}{}
     \item{\link[postNet]{ptn_geneList}}{}
     \item{\link[postNet]{ptn_GO}}{}
     \item{\link[postNet]{ptn_GSEA}}{}
     \item{\link[postNet]{ptn_id}}{}
     \item{\link[postNet]{ptn_miRNA_analysis}}{}
     \item{\link[postNet]{ptn_miRNA_to_gene}}{}
     \item{\link[postNet]{ptn_model}}{}
     \item{\link[postNet]{ptn_motifgeneList}}{}
     \item{\link[postNet]{ptn_motifSelection}}{}
     \item{\link[postNet]{ptn_selectedFeatures}}{}
     \item{\link[postNet]{ptn_selection}}{}
     \item{\link[postNet]{ptn_sequences}}{}
     \item{\link[postNet]{ptn_species}}{}
     \item{\link[postNet]{ptn_version}}{}
     }
}

\examples{

# -----------------------------------------------------------------------
# Example 1: Using an anota2seq object with built-in RefSeq annotations
# -----------------------------------------------------------------------

# load example data with an Anota2seqDataSet object:
data("postNetExample")

# This postNet analysis will allow identification of mRNA features, 
# and comparison between translationally activated and suppressed genes.

# Feature integration analysis will model changes in translation
# efficiency.

# Initialize the postNetData object, where "ads" is an Anota2seqDataSet object
# resulting from an anota2seqRun call:
    ptn <- postNetStart(
      ads = ads,
      regulation = c("translationUp","translationDown"),
      contrast = c(1,1),
      regulationGen = "translation",
      contrastSel = 1,
      selection = "random",
      setSeed = 123,   # ensures reproducibility of random isoform selection
      source = "load",
      species = "human"
    )

# -----------------------------------------------------------------------
# Example 2: Using custom gene lists and background
# -----------------------------------------------------------------------

# Instead of using the results of an anota2seq analysis, gene sets of interest
# and custom background can also be provided, along with a regulation effect measurement. 

# Genes of interest should be provided in a named list.
 myGenes <- ptn_getExampleGeneList(ads)
 str(myGenes)
 
# All gene IDs in the list should be present in the background.
 myBg <- ptn_getExampleBg(ads)
 str(myBg)
 
# The regulation effect measurement must be named with the same gene IDs
# present in the background (or the gene list if no custom background is provided).
 myEffect <- ptn_getExampleEffect(ads)
 str(myEffect)

# Initialize the postNetData object with custom gene lists:
 ptn <- postNetStart(
   geneList = myGenes,
   geneListcolours = c("#7FB7BE","#DB7F67"),
   customBg = myBg,
   effectMeasure = myEffect,
   selection = "random",
   setSeed = 123,
   source = "load", 
   species = "human"
 )

\dontrun{
# -----------------------------------------------------------------------
# Example 3: Creating annotation from local RefSeq files
# -----------------------------------------------------------------------

#  Example RefSeq annotation files for mouse downloaded from the NCBI database:
#  - "myMouse_rna.gbff.gz"
#  - "myMouse_rna.fa.gz"
#  - "myMouse_genomic.gff.gz"
 
# Initialize the postNetData object with custom gene lists and RefSeq annotations,
# and considering the longest mRNA isoform for genes with multiple isoforms:
 ptn <- postNetStart(
   geneList = myGenes,
   geneListcolours = c("#F2A104","#00743F"),
   customBg = myBg,
   effectMeasure = myEffect,
   selection = "longest",
   source = "createFromSourceFiles",
   species = "mouse",
   rna_gbff_file = "myMouse_rna.gbff.gz",
   rna_fa_file = "myMouse_rna.fa.gz",
   genomic_gff_file = "myMouse_genomic.gff.gz"
 )
}

# -----------------------------------------------------------------------
# Example 4: Retrieving data from the postNetData object
# -----------------------------------------------------------------------

# Inputs and results stored in the slots of the postNetData object can be
# easily retrieved with the use of ptn helper functions.

# To extract reference annotation sequences for 5'UTRs from a postNetData object ptn:
 myUTR5seqs <- ptn_sequences(
    x = ptn,
    region = "UTR5"
 )

 str(myUTR5seqs)
 
# To extract the RefSeq annotation release version from postNetData object ptn (if using source = "load"):
 myVersion <- ptn_version(
    x = ptn)

 myVersion

}

\seealso{
\code{\link[anota2seq]{anota2seqRun}} for obtaining the anota2seq \code{ads} object. \cr
\code{\link[postNet]{checkAvailableVersions}} to see the version numbers of the in-build reference sequence annotations. \cr
\code{\link[postNet]{featureIntegration}} for requirements of postNet feature integration analysis.
}

\references{
\url{https://www.ncbi.nlm.nih.gov/refseq/} for RefSeq documentation.  
\url{https://bioconductor.org/packages/release/bioc/vignettes/anota2seq/inst/doc/anota2seq.pdf} for \pkg{anota2seq} details.
}
