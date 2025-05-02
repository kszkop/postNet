\name{plotSignatures}
\alias{plotSignatures}
\title{Analysis of regulation of gene signatures in a \code{postNetData} object
}

\description{
}

\usage{
plotSignatures(ptn, 
              signatureList, 
              signature_colours = NULL, 
              xlim = c(-2, 2), 
              tableCex = 1, 
              pdfName = NULL)
}

\arguments{
  \item{ptn}{A \code{postNetData} object
  }
  \item{signatureList}{A named list of vectors containing gene IDs for the signatures of interest to be examined. Note that several
  signatures of translational regulation are provided with the package for both human and mouse. These can be retrieved using the
  \link[postNet]{retrieveSignatures} function.
  }
  \item{signature_colours}{An optional character vector specifying colours for plotting each signature. Note that \code{signature_colours} must be the same length as
  \code{signatureList}. The default is \code{NULL}. 
  }
  \item{xlim}{A numeric vector specifying the x-axis limits for the ecdf plots. Default is \code{c(-2,2)}.
  }
  \item{tableCex}{
  }
  \item{pdfName}{
  }
}

\details{
}

\value{No value is returned. Graphical outputs are generated as PDF files. 
}

\seealso{
\link[postNet]{signaturesHeatmap} \cr
\link[postNet]{plotSignatures_ads} \cr
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
 
# load signature data:
 data('humanSignatures')
 
# Examine the regulation of mTOR-sensitive transcripts in the results of the anota2seq analysis 
plotSignatures(ptn = ptn,
                               signatureList = humanSignatures[c('Gandin_etal_2016_mTOR_transUp','Gandin_etal_2016_mTOR_transDown')],
                               signature_colours = c('red','blue'), 
                               xlim = c(-3,3), 
                               tableCex = 1,
                               pdfName = 'mTORsignatures')

}
