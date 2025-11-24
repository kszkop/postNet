\name{lengthAnalysis}
\alias{lengthAnalysis}
\title{Caculate the length of mRNA sequence regions and compare between gene sets
}

\description{The \code{lengthAnalysis} function calculates the log2-transformed length of different regions of mRNA 
sequences (5'UTR, 3'UTR, and CDS) for each gene. In addition, the length of different sequence regions for gene sets 
of interest can be compared against background, or other genes sets. Several options are available for plotting the 
comparisons. 
}

\usage{
lengthAnalysis(ptn, 
              region, 
              comparisons = NULL,
              plotOut = TRUE, 
              plotType = "boxplot", 
              pdfName = NULL)
}

\arguments{
  \item{ptn}{A \code{postNetData} object.
  }
  \item{region}{A character vector specifying the sequence region(s) to be analyzed. This can be \code{"UTR5"}, \code{"UTR3"}, \code{"CDS"}, or any combination of the regions.
  }
  \item{comparisons}{A list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the 
  log2-transformed lengths for the sequence regions of interest will be returned, and statistical 
  comparisons between gene sets will not be performed. The default is \code{TRUE}.
  }
  \item{plotType}{If \code{plotOut} is \code{TRUE}, indicate the type of output plot to generate. Options are,
  \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"boxplot"}. 
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{A two-sided Wilcoxon Rank Sum test is performed to identify significant differences in mRNA sequence region length between gene sets of interest, or against the background gene set.
}

\value{A named list of vectors for selected mRNA sequence regions with the calculated log2-transformed lengths for each gene. The elements in each vector are named with the gene ID. The \code{lengthAnalysis} function can also return PDF files of output plots with statistical comparisons between gene sets of interest.
}

\examples{
# load example data:
 ptn <- postNetExample$ptn

# Calculate the length of the 5'UTR, 3'UTR and coding sequence for each gene, and compare lengths between
# translationUp genes vs. background, translationDown genes vs. background, and translationUp vs. translationDown genes
 len <- lengthAnalysis(ptn = ptn, 
                      region = c("UTR5","CDS","UTR3"), 
                      comparisons = list(c(0,1),c(0,2),c(1,2)),
                      plotOut = TRUE,
                      plotType = 'boxplot',
                      pdfName = "example")
                      
 str(len)

}
