\docType{methods}
\name{ptn_codonAnalysis}
\alias{ptn_codonAnalysis,postNetData-method}
\title{
A Capitalized Title (ideally limited to 65 characters)
}

\description{
}

\usage{
\S4method{ptn_codonAnalysis}{postNetData}(x)
}

\arguments{
  \item{x}{
}
}

\value{

}

\seealso{
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