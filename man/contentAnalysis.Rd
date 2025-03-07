\name{contentAnalysis}
\alias{contentAnalysis}
\title{Caculate and compare the nucleotide content of mRNA sequence regions between gene sets 
}

\description{Nucleotide content of desired region for each gene is calculated. If plotOut is TRUE, all indicated comparisons between different sets of genes are performed and plotted. 
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
  \item{subregion}{Optionally it is possible to select a subregion by indicating a number of nucleotides 
  from the start if positive or end if negative. 
  }
  \item{subregionSel}{It is required if the subregion is provided. Choose whether the indicated subregion
  should be "select" or "exclude".
  }
  \item{comparisons}{A list of numeric vectors indicating which gene sets in \code{regulation} or 
  \code{geneList} should be used in pair-wise comparisons. The numeric vectors should correspond to the
  indexes of the gene sets, where 0 denotes the background. For example, \code{list(c(0,2),c(1,2))} 
  indicates that the second gene set in \code{regulation} or \code{geneList} will be compared to the
  background (0), and the first and second gene sets will be compared to each other. The default is
  \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the 
  log2-transformed lengths for the sequence regions of interest will be returned, and comparisons between
  gene sets will not be performed. The default is \code{TRUE}. 
  }
  \item{plotType}{If \code{plotOut} is \code{TRUE}, indicate the type of output plot to generate. Options 
  are, \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"boxplot"}. 
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{A two-sided Wilcoxon Rank Sum test is performed to identify significant differences in mRNA sequence region length between gene sets of interest, or against the background gene set.
}

\value{List of vectors for selected regions with calculated desired nucleotide content for each gene. The function will also output plots with desired comparisons.
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
 

# Calculate the GC of the 5'UTR, 3'UTR and coding sequence for each gene, and compare between
# translationUp genes vs. background, translationDown genes vs. background, and translationUp vs. translationDown genes
 content <- lengthAnalysis(ptn = ptn, 
                      region = c("UTR5","CDS","UTR3"), 
                      comparisons = list(c(0,1),c(0,2),c(1,2)),
                      plotOut = TRUE,
                      plotType = 'boxplot',
                      pdfName = "example")
                      
 str(content)
}
