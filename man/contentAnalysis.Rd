\name{contentAnalysis}
\alias{contentAnalysis}
\title{Caculate the nucleotide content of mRNA sequence regions and compare between gene sets 
}

\description{The \code{contentAnalysis} function calculates the percentage of nucleotide content of mRNA sequences (5'UTR, 3'UTR, and CDS) for each gene. In addition, the nucleotide content of different sequence regions for gene sets of interest can be compared against background, or other genes sets. Several options are available for plotting the comparisons.
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
  \item{ptn}{a \code{postNetData} object.
  }
  \item{contentIn}{A vector with desired content combinations of A,T,G,C
  }
  \item{region}{A character vector specifying the sequence region(s) of interest. This can be \code{"UTR5"},
  \code{"UTR3"}, \code{"CDS"}, or any combination of the regions.
  }
  \item{subregion}{Optionally, it is possible to specify a more specific subregion of the sequences to either
  select or exclude from the analysis. This is done by providing a numeric value indicating the number of
  nucleotides from either the start of the sequence region if positive, or the end if negative. For example,
  \code{subregion = 50}, would limite the analysis to the first 50 nucleotides of the sequence region(s)
  specified by \code{region}, while \code{subregion = -50}, would remove the last 50 nucleotides of the
  sequence region(s). The default is \code{NULL}.  
  }
  \item{subregionSel}{If a value is provided for \code{subregion}, indicate whether the subregion should be
  considered or excluded by specifying either \code{"select"} or \code{"exclude"}. The default is \code{NULL}.
  }
  \item{comparisons}{A list of numeric vectors indicating which gene sets in \code{regulation} or 
  \code{geneList} should be used in pair-wise comparisons. The numeric vectors should correspond to the
  indexes of the gene sets, where 0 denotes the background. For example, \code{list(c(0,2),c(1,2))} 
  indicates that the second gene set in \code{regulation} or \code{geneList} will be compared to the
  background (0), and the first and second gene sets will be compared to each other. The default is
  \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the 
  nucleotide content for the sequence regions of interest will be returned, and statistical 
  comparisons between gene sets will not be performed. The default is \code{TRUE}. 
  }
  \item{plotType}{If \code{plotOut} is \code{TRUE}, indicate the type of output plot to generate. Options 
  are, \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"boxplot"}. 
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{A two-sided Wilcoxon Rank Sum test is performed to identify significant differences in mRNA sequence region nucleotide content between gene sets of interest, or against the background gene set.
}

\value{A named list of vectors for selected mRNA sequence regions with the calculated nucleotide content for each gene. The elements in each vector are named with the gene ID. The \code{contentAnalysis} function can also return PDF files of output plots with statistical comparisons between gene sets of interest.
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
 

# Calculate the GC of the first 50 nucleotidfes of the 5'UTR for each gene, and compare between
# translationUp genes vs. background, translationDown genes vs. background, and translationUp vs. translationDown genes
 content_UTR5_first50 <- contentAnalysis(ptn = ptn, 
                           region = c("UTR5"),
                           subregion = 50,
                           subregionSel = "select",
                           comparisons = list(c(0,1),c(0,2),c(1,2)), 
                           contentIn = c('GC'),
                           plotOut = TRUE,
                           plotType = 'ecdf',
                           pdfName = 'first50')
                      
 str(content_UTR5_first50)
}
