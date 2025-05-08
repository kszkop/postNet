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
  \item{tableCex}{A numeric value specifying the size of the text in the table of statistical results that appears at the top of the ecdf plots.
  Default is \code{1}.
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{Gene signatures of interest are visualized as using empirical cumulative distribution functions (ECDFs) of log2 fold changes in the regulatory effect
measurement in a \code{postNetData} object. ECDFs for genes belonging to the signatures of interest are compare to those for all other genes 
(background). Differences in the regulatory effect measurement distributions are calculated at the quantiles. Significant directional shifts in the ECDFs for gene 
signature versus background are identified using a Wilcoxon rank-sum test.
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
