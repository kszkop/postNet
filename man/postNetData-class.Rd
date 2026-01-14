\name{postNetData-class}
\docType{class}
\alias{postNetData-class}

\title{S4 Class \code{postNetData}: Core data container for postNet analyses
}

\description{The \code{postNetData} class stores the core data structure used throughout the postNet analysis. It contains input gene lists, sequence annotations and source information, and results from codon, motif, GSEA, GAGE, GO term, and feature integration analyses. This object is created by the \code{\link[postNet]{postNetStart}} function and serves as the input/output container for most downstream functions.
}

\section{Objects from the Class}{
Objects of this class are typically created using the constructor function \code{\link[postNet]{postNetStart}}. Alternatively, they can be created manually using \code{new("postNetData", ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{species}:}{Object of class \code{"characterOrNULL"} specifying the species associated with the sequence annotations and input gene lists, e.g. \code{"human"}, \code{"mouse"}. Used to load species-specific reference sequences.
 }
    \item{\code{version}:}{Object of class \code{"characterOrNULL"} storing the version identifier for the RefSeq sequence annotations loaded.
}
    \item{\code{selection}:}{Object of class \code{"character"} indicating the method used to select representative transcript isoforms, e.g. \code{"random"}, \code{"shortest"}, or \code{"longest"}.
}
    \item{\code{seed}:}{Object of class \code{"ANY"} storing the integer used to set the seed for random sampling if \code{selection = "random"}
}
    \item{\code{annot}:}{Object of class \code{"postNetAnnot"} storing the annotated sequences for each mRNA transcript region, \code{"UTR5"}, \code{"CDS"}, \code{"CCDS"}, and \code{"UTR3"}.
}
    \item{\code{dataIn}:}{Object of class \code{"postNetDataIn"} containing the user-defined input gene lists of interest, background gene set, regulatory effect measurement, and selected colours for plotting.
}
    \item{\code{features}:}{Object of class \code{"dataframeOrNULL"} with a \code{data.frame} containing the features for each gene (e.g., UTR lengths, GC content, codon or motif counts) that were used as input for stepwise regression and/or random forest models during \code{\link[postNet]{featureIntegration}} analysis.
}
    \item{\code{analysis}:}{Object of class \code{"postNetAnalysis"} storing the results from analysis functions including \code{codonUsage}, \code{motifAnalysis}, \code{gseaAnalysis}, \code{gageAnalysis}, \code{goAnalysis}, \code{miRNAanalysis} and \code{featureIntegration}. Structured into nested S4 objects.
}
  }
}
\section{Methods}{The following S4 accessor and utility methods are available for \code{postNetData} objects:
  \describe{
    \item{ptn_background}{\code{signature(x = "postNetData")}: Retrieve the background gene set.
    }
    \item{ptn_codonAnalysis}{\code{signature(x = "postNetData")}: Access results from the codon or amino acid usage analysis.}
    \item{ptn_codonSelection}{\code{signature(x = "postNetData")}: Retrieve the selected significantly enriched/depleted codons or amino acids.}
    \item{ptn_colours}{\code{signature(x = "postNetData")}: Return the colours that will be used for visualizations.}
    \item{ptn_dataIn}{\code{signature(x = "postNetData")}: Access the user-defined background gene set, gene lists of interest, regulatory effect measurements, and colours for visualizations.}
    \item{ptn_effect}{\code{signature(x = "postNetData")}: Retrieve the user-supplied regulatory effect measurement (e.g., Log2 Fold Changes).}
    \item{ptn_features}{\code{signature(x = "postNetData")}: Extract the \code{features} \code{data.frame} used as input for stepwise regression and random forest analyses.}
    \item{ptn_geneID}{\code{signature(x = "postNetData")}: Return a character vector of gene IDs for the genes included in the gene sets of interest and background gene sets (all genes used in the analysis).}
    \item{ptn_geneList}{\code{signature(x = "postNetData")}: Access the list of user-supplied gene sets of interest used in analyses.}
    \item{ptn_id}{\code{signature(x = "postNetData")}: Return a character vector of transcript IDs for the transcripts included in the gene sets of interest and background gene sets (all transcripts used in the analysis).}
    \item{ptn_motifgeneList}{\code{signature(x = "postNetData")}: Access the results and additional information on the de novo identified motifs for the specified gene list in \link[universalmotif]{universalmotif} format.}
    \item{ptn_motifSelection}{\code{signature(x = "postNetData")}: Return a character vector of de novo identified motifs in the gene sets of interest that have enrichment p-values below a user-specified threshold.}
    \item{ptn_networkGraph}{\code{signature(x = "postNetData")}: Retrieve the network structure in igraph format generated during linear stepwise regression modelling with the \link[postNet]{featureIntegration} function.}
    \item{ptn_selection}{\code{signature(x = "postNetData")}: Return the user-defined transcript isoform selection strategy (e.g., \code{random}, \code{shortest}, \code{longest}).}
    \item{ptn_sequences}{\code{signature(x = "postNetData")}: Retrieve a character vector containing the reference annotation sequences corresponding to the specified transcript region.}
    \item{ptn_species}{\code{signature(x = "postNetData")}: Return the species associated with the analysis.}
    \item{ptn_version}{\code{signature(x = "postNetData")}: Return the release version of the RefSeq annotations used to generate the reference sequences.}
     \item{ptn_GAGE}{\code{signature(x = "postNetData")}: Access the results of GAGE analysis.}
     \item{ptn_GSEA}{\code{signature(x = "postNetData")}: Access the results of GSEA analysis.}
     \item{ptn_GO}{\code{signature(x = "postNetData")}: Access the results of GO term analysis.}
      \item{ptn_miRNA_analysis}{\code{signature(x = "postNetData")}: Access the results of GAGE-based miRNA enrichment analysis.}
       \item{ptn_miRNA_to_gene}{\code{signature(x = "postNetData")}: Return a list of genes that are predicted to be targeted by the miRNA that passed the filtering thresholds in miRNA analysis.}
	 }
}

\seealso{
\link[postNet]{postNetStart} \cr
\link[postNet]{motifAnalysis} \cr
\link[postNet]{codonUsage} \cr
\link[postNet]{miRNAanalysis} \cr
\link[postNet]{gseaAnalysis} \cr
\link[postNet]{gageAnalysis} \cr
\link[postNet]{goAnalysis} \cr
\link[postNet]{featureIntegration} \cr
\link[postNet]{plotFeaturesMap} \cr
\link[postNet]{ptn_background} \cr
\link[postNet]{ptn_check_models} \cr
\link[postNet]{ptn_codonAnalysis} \cr
\link[postNet]{ptn_codonSelection} \cr
\link[postNet]{ptn_colours} \cr
\link[postNet]{ptn_dataIn} \cr
\link[postNet]{ptn_effect} \cr
\link[postNet]{ptn_features} \cr
\link[postNet]{ptn_GAGE} \cr
\link[postNet]{ptn_geneID} \cr
\link[postNet]{ptn_geneList} \cr
\link[postNet]{ptn_GO} \cr
\link[postNet]{ptn_GSEA} \cr
\link[postNet]{ptn_id} \cr
\link[postNet]{ptn_miRNA_analysis} \cr
\link[postNet]{ptn_miRNA_to_gene} \cr
\link[postNet]{ptn_model} \cr
\link[postNet]{ptn_motifgeneList} \cr
\link[postNet]{ptn_motifSelection} \cr
\link[postNet]{ptn_networkGraph} \cr
\link[postNet]{ptn_selectedFeatures} \cr
\link[postNet]{ptn_selection} \cr
\link[postNet]{ptn_sequences} \cr
\link[postNet]{ptn_species} \cr
\link[postNet]{ptn_version}
}

\examples{
showClass("postNetData")
}

\keyword{classes}
