\docType{methods}
\name{ptn_codonAnalysis}
\alias{ptn_codonAnalysis,postNetData-method}
\title{Accessor for the "codonAnalysis" slot of a \code{postNetData} object
}

\description{Retrieve the "codonAnalysis" slot from a \code{postNetData} object, which holds the results of of the codon usage analysis produced using the \link[postNet]{codonUsage} function.
}

\usage{
\S4method{ptn_codonAnalysis}{postNetData}(x)
}

\arguments{
  \item{x}{A \code{postNetData} object.
}
}

\value{The value returned by the \code{ptn_codonAnalysis} function is a data.frame with the following columns: geneID, codon, AA, count, frequency, AACountPerGene, and relative_frequency. Here, frequency corresponds to the number of codons relative to all codons in the gene, and relative_frequency corresponds to the number of codons relative to all synonymous codons in the gene.
}

\seealso{
\link[postNet]{codonUsage} \cr
\link[postNet]{ptn_codonSelection}
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
 
# Run codon usage analysis of single codons
 ptn <- codonUsage(ptn = ptn, 
          annotType = 'ptnCDS', 
          sourceSeq = "load", 
          analysis = 'codon', 
          codonN = 1, 
          pAdj = 0.01, 
          rem5 = TRUE, 
          plotHeatmap = FALSE, 
          thresOddsUp = 0.3, 
          thresFreqUp = 0.3, 
          thresOddsDown = 0.3, 
          thresFreqDown = 0.3, 
          subregion = NULL, 
          subregionSel = NULL, 
          comparisons = list(c(1,2)), 
          plotType_index = "violin", 
          pdfName = 'codons')
          
# Get the results of the codon usage analysis
 
 codonResults <- ptn_codonAnalysis(ptn)
 str(codonResults)
}
