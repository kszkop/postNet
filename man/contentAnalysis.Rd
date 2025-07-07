\name{contentAnalysis}
\alias{contentAnalysis}
\title{Calculate nucleotide content of mRNA sequence regions and compare between gene sets 
}

\description{The \code{contentAnalysis} function calculates the nucleotide content (as a percentage) of mRNA sequences (5'UTR, 3'UTR, and CDS) for each gene. Nucleotide content of different sequence regions can be compared between gene sets of interest, or against background. Several options are available for plotting the comparisons.
}

\usage{
contentAnalysis(ptn,
                contentIn, 
                region, 
                subregion = NULL, 
                subregionSel = NULL, 
                comparisons = NULL, 
                plotOut = TRUE, 
                plotType = "boxplot", 
                pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }
  \item{contentIn}{A character vector specifying one or more nucleotide(s) or nucleotide combinations to quantify. These can be any selection or combination of A, T, G, or C. For example \code{c("GC")}, or \code{c("G", "A")}. When \code{region} includes \code{"CDS"}, nucleotide content at specific codon positions can also be quantified. For example, \code{c("GC3")} calculates GC content in the third codon position, and \code{c("A12")} calculates the percentage of A content in the first and second codon positions. 
  }
  \item{region}{A character vector specifying the sequence region(s) to be analyzed. This can be \code{"UTR5"}, \code{"UTR3"}, \code{"CDS"}, or any combination of the regions.
  }
  \item{subregion}{Optionally, specify a more specific subregion of the sequences to either select or exclude from the analysis. Provide an integer indicating the number of nucleotides to include or exclude from either the start of the sequence region if positive, or the end if negative. For example, \code{subregion = 50}, denotes the first 50 nucleotides of the sequence region(s) specified by \code{region}, while \code{subregion = -50}, denotes the last 50 nucleotides of the sequence region(s). The default is \code{NULL}.  
  }
  \item{subregionSel}{If a value is provided for \code{subregion}, indicate whether the subregion should be selected (limiting the analysis to this region) or excluded by specifying either \code{"select"} or \code{"exclude"}. The default is \code{NULL}.
  }
  \item{comparisons}{A list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the nucleotide content for the sequence regions of interest will be returned, and statistical comparisons between gene sets will not be performed. The default is \code{TRUE}. 
  }
  \item{plotType}{If \code{plotOut} is \code{TRUE}, indicate the type of output plot to generate. Options are \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"boxplot"}. 
  }
  \item{pdfName}{Name to be appended to output PDF files. The default is \code{NULL}.  
  }
}

\details{When \code{plotOut = TRUE}, a two-sided Wilcoxon Rank Sum test is performed to identify significant differences in mRNA sequence region nucleotide content between gene sets of interest, and/or against the background gene set.
}

\value{A named list of vectors for selected mRNA sequence regions with the calculated nucleotide content for each gene. This list can be used as input for the downstream \link[postNet]{featureIntegration} analysis. The elements in each vector are named with the gene ID. The \code{contentAnalysis} function can also return PDF files of output plots with statistical comparisons between gene sets of interest.
}

\examples{
# load example data:
 data(postNetExample)
    
# Initialize the postNetData object with custom gene lists:
 myGenes <- postNetExample$geneList
 myBg <- postNetExample$background
 myEffect <- postNetExample$effect

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
 

# Calculate the GC of the first 50 nucleotides of the 5'UTR for each gene, and compare between
# translationUp genes vs. background, translationDown genes vs. background, and translationUp vs. translationDown genes
 content_UTR5_first50 <- contentAnalysis(ptn = ptn, 
                           region = c("UTR5"),
                           subregion = 50,
                           subregionSel = "select",
                           comparisons = list(c(0,1),c(0,2),c(1,2)), 
                           contentIn = c("GC"),
                           plotOut = TRUE,
                           plotType = "ecdf",
                           pdfName = "first50")
                      
 str(content_UTR5_first50)
}
