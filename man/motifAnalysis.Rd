\name{motifAnalysis}
\alias{motifAnalysis}
\title{De-novo motif discovery using STREME and a \code{postNetData} object
}

\description{The \code{motifAnalysis} function applies the STREME method from the MEME-Suite to discover de-novo ungapped motifs (recurring, fixed-length patterns) that are enriched in the selected sequence region(s) of gene sets of interest relative to control sequences (background). 
}

\usage{
motifAnalysis(ptn,
              stremeThreshold = 0.05, 
              minwidth = 6, 
              memePath = NULL, 
              seqType = "dna", 
              region, 
              subregion = NULL, 
              subregionSel = NULL)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
  \item{stremeThreshold}{A numeric value specifying the enrichment p-value threshold for motif selection. 
  The default is 0.05. 
  }
  \item{minwidth}{A numeric value specifying the minimal width (in nucleotides) of the motif. The default is 6.
  }
  \item{memePath}{The path to the STREME executables in "meme/bin". Note that the MEME-Suite must be installed. 
  }
  \item{seqType}{A string speicfying the type of sequence being provided to search for motifs. The options are,
  \code{"dna"}, \code{"rna"}, or \code{"protein"}. The default is \code{"dna"}.
  }
  \item{region}{A character vector specifying the sequence region(s) of interest. This can be \code{"UTR5"},
  \code{"UTR3"}, \code{"CDS"}, or any combination of the regions.
  }
  \item{subregion}{Optionally, it is possible to specify a more specific subregion of the sequences to either
  select or exclude from the analysis. This is done by providing a numeric value indicating the number of
  nucleotides from either the start of the sequence region if positive, or the end if negative. For example,
  \code{subregion = 50}, would denote the first 50 nucleotides of the sequence region(s) specified by \code{region},
  while \code{subregion = -50}, would denote the last 50 nucleotides of the sequence region(s). The
  default is \code{NULL}. 
  }
  \item{subregionSel}{If a value is provided for \code{subregion}, indicate whether the subregion should be
  selected (limiting the analysis to this region) or excluded by specifying either \code{"select"} or 
  \code{"exclude"}. The default is \code{NULL}.
  }
}

\value{The results returned by the \code{motifAnalysis} function are compiled into an S4 object of class \code{postNetMotifs} and stored in the "motifs" slot of the \code{postNetData} object. Within the motif slot, results for each sequence region are provided as lists. The "motifSelection" list element corresponds to identified motifs in the gene lists of interest that have enrichment p-values below the threshold set with the \code{stremeThreshold} parameter. Additional list elements are provided for each gene set of interested with the selected results with further information about the motifs in \link[universalmotif]{universalmotif} format. 
}

\details{The \code{motifAnalysis} function applies the \link[memes]{runStreme} function from the \link[memes]{memes} package. However, the MEME Suite software tools must also be installed locally, and the path to the executables must be provided. The MEME Suite can be downloaded here: \url{https://meme-suite.org/meme/doc/download.html}. Installation instructions can be found here: \url{https://meme-suite.org/meme/doc/install.html?man_type=web}. 

Motifs identified using \code{motifAnalysis} can be quantified and included in downstream \link[postNet]{featureIntegration} using the \link[postNet]{contentMotifs} function. Significantly enriched motifs can be extracted using the S4 method \link[postNet]{ptn_motifSelection} function, and \link[universalmotif]{universalmotif} results for each gene set included in \code{geneList} can be retrieved using the \link[postNet]{ptn_motifgeneList} function.   

If you use the \code{motifAnalysis} function in your analysis, please appropriately cite The MEME Suite and the STREME tool publications provided below. Licensing The MEME Suite is free for non-profit use. However, for-profit users should contact OIC-MEMESuite@ucsd.edu. See \url{http://meme-suite.org/doc/copyright.html} for details}

\references{

Timothy L. Bailey, James Johnson, Charles E. Grant, William S. Noble, "The MEME Suite", Nucleic Acids Research, 43(W1):W39-W49, 2015.

Timothy L. Bailey, "STREME: Accurate and versatile sequence motif discovery", Bioinformatics, 2021. https://doi.org/10.1093/bioinformatics/btab203
}

\seealso{
\link[memes]{memes} \cr
\link[memes]{runStreme} \cr
\link[universalmotif]{universalmotif} \cr
\link[postNet]{contentMotifs}
\link[postNet]{ptn_motifSelection}
\link[postNet]{ptn_motifGeneList}
}

\examples{
# Note that as users must provide the path to the MEME Suite executables, it is not possile to
# create a fully executable example as this path will vary depending on many factors. An example of
# how to run this function is provided below, and can be updated with the correct memePath argument. 

\dontrun{
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
 
 ptn <- motifAnalysis(ptn = ptn, 
                     stremeThreshold = 0.05,
                     minwidth = 6,
                     memePath = "/meme/bin", 
                     region = c('UTR5','CDS','UTR3'),
                     subregion = NULL,
                     subregionSel = NULL)
                     
# Extract the significantly enriched motifs:
 myMotifs <- ptn_motifSelection(ptn = ptn,
                                region = 'UTR5')
                                
 str(myMotifs)

# Extract full motif results from one gene set of interest:                                
  motifsTransUp <- ptn_motifSelection(ptn = ptn,
                                      region = 'UTR5',
                                      geneList = 'translationUp')
                                
 str(motifsTransUp)
}

}

