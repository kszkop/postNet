\name{plotSignatures_ads}
\alias{plotSignatures_ads}
\title{Analysis of regulation of gene signatures in an \code{Anota2seqDataSet} object
}

\description{The \code{plotSignatures_ads} function assesses changes in the regulation of known gene signatures of interest in an 
\code{Anota2seqDataSet} object from an \pkg{anota2seq} analysis. Outputs include statistical analyses and visualizations of signature
regulation. 
}

\usage{
plotSignatures_ads(ads, 
                  contrast, 
                  dataName, 
                  effects_names = c("Total mRNA Log2FC", "Polysome-associated mRNA Log2FC", "Buffering Log2FC", "Translation Log2FC"), 
                  signatureList, 
                  generalName, 
                  signature_colours, 
                  xlim = c(-2, 2), 
                  scatterXY = NULL, 
                  tableCex = 1, 
                  pdfName = NULL)
}

\arguments{
  \item{ads}{An S4 object of class \code{Anota2seqDataSet} (resulting from the \link[anota2seq]{anota2seqRun} function in the \pkg{anota2seq} package). 
  }
  \item{contrast}{An integer selecting the \pkg{anota2seq} contrast to retrieve regulatory effect measurements from. 
  For more details on contrasts, please see the \pkg{anota2seq} vignette. 
  }
  \item{dataName}{A string specifying the name that will be given to selected contrast from the \code{Anota2seqDataSet} object. For example: \code{"Test vs. DMSO"}. 
  }
  \item{effects_names}{A character vector specifying axis names for the fold change scatterplot and ecdf plots. The default is, 
  \code{c("Total mRNA Log2FC", "Polysome-associated mRNA Log2FC", "Buffering Log2FC", "Translation Log2FC")}. 
  }
  \item{signatureList}{A named list of vectors containing gene IDs for the signatures of interest to be examined. Note that several
  signatures of translational regulation are provided with the package for both human and mouse (see details and example below). 
  }
  \item{generalName}{A character vector specifying the names of the signatures to be plotted. 
  }
  \item{signature_colours}{A character vector specifying colours for plotting each signature. Note that \code{signature_colours} must be the same length as
  \code{signatureList}. 
  }
  \item{xlim}{A numeric vector specifying the x-axis limits for the ecdf plots. Default is \code{c(-2,2)}.
  }
  \item{scatterXY}{A numeric value specifying the scale of the x and y axes for the fold change scatterplot. For example, \code{scatterXY = 3}, would indicate that both
  the x and y axis will have a Log2FC range from -3 to 3. The default is \code{NULL}. 
  }
  \item{tableCex}{A numeric value specifying the size of the text in the table of statistical results that appears at the top of the ecdf plots.
  Default is \code{1}. 
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{Gene signatures of interest are visualized in a scatterplot of log2 fold changes in polysome-associated and total mRNA. The empirical 
cumulative distribution functions (ECDFs) of log2 fold changes in polysome-associated mRNA, total mRNA, translation, and buffering (offsetting) 
are also plotted independently. Fold change ECDFs for genes belonging to the signatures of interest are compared to those for all other genes 
(background). Differences in the fold change distributions are calculated at the quantiles. Significant directional shifts in the ECDFs for gene 
signature versus background are identified using a Wilcoxon rank-sum test.
}

\value{No value is returned. Graphical outputs are generated as PDF files. 
}

\seealso{
\link[anota2seq]{anota2seqRun} \cr
\link[postNet]{signaturesHeatmap} \cr
\link[postNet]{plotSignatures} \cr
}

\examples{
# load example data:
 data("postNetExample", package = "postNetParcel")

# Initialize Anota2seqDataSet (see anota2seq vignette for details)
 ads <- anota2seq::anota2seqDataSetFromMatrix(
    dataP = postNetExample$ads_data$dataP,
    dataT = postNetExample$ads_data$dataT,
    phenoVec = postNetExample$ads_data$phenoVec,
    batchVec = c(1, 2, 3, 4, 1, 2, 3, 4),
    dataType = "RNAseq",
    normalize = FALSE)
    
# Run an anota2seq analysis:
# Note that the quality control and residual outlier testing are not 
# performed to limit the running time of this example. For full details 
# on running an analysis please see the anota2seq vignette and help manual.
 ads <- anota2seq::anota2seqRun(ads,
    performQC = FALSE, 
    performROT = FALSE, 
    useProgBar = FALSE)
    
# load signature data:
 humanSignatures <- get_signatures("human")
 
# Examine the regulation of mTOR-sensitive transcripts in the results of the anota2seq analysis 
plotSignatures_ads(ads = ads,
                               contrast = 1,
                               dataName = "Osmosis_1h",
                               signatureList = humanSignatures[c('Gandin_etal_2016_mTOR_transUp',
                                                                'Gandin_etal_2016_mTOR_transDown')],
                               generalName = "mTOR-sensitive translation",
                               signature_colours = c('red','blue'), 
                               xlim = c(-3,3), 
                               scatterXY = 4,
                               tableCex = 0.65,
                               pdfName = 'mTORsignatures')

}
