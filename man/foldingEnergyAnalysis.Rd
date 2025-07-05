\name{foldingEnergyAnalysis}
\alias{foldingEnergyAnalysis}
\title{Compare mRNA folding energy between gene sets  
}

\description{The \code{foldingEnergyAnalysis} function compares folding energies for each mRNA sequence region. Pre-calculated folding energies for reference sequences are available for certain RefSeq releases (for human and mouse). These values were calculated using the mfold algorithm. Alternatively, user-supplied values can be used in analyses. In addition, folding energies for gene sets of interest can be compared against background, or other gene sets. 
Several options are available for plotting the comparisons. 
}

\usage{
foldingEnergyAnalysis(ptn, 
                      sourceFE = "load", 
                      customFileFE = NULL, 
                      residFE = FALSE, 
                      region, 
                      comparisons = NULL, 
                      plotOut = TRUE, 
                      plotType = "ecdf", 
                      pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
  }
  \item{sourceFE}{A string specifying the source of the folding energy values. The options are \code{"load"}, or \code{"custom"}. If \code{"load"}, pre-calculated folding energies for the selected species ("human" or "mouse") will be used. If \code{"custom"}, a custom file with pre-calculated folding energies must be provided using the \code{customFileFE} argument. The default is \code{"load"}. 
  }
  \item{customFileFE}{If \code{sourceFE = "custom"}, provide a string specifying the name and file path to the custom folding energy file. The custom file with calculated folding energies must be provided in a TAB delimited format with three columns: transcript ID, folding energy, and length of the sequence. Note that folding energies for different sequence regions (e.g. 5'UTR and 3'UTR) must be provided separately.
  }
  \item{residFE}{Logical indicating if the folding energies should be corrected for the length of the sequences. If \code{TRUE}, the values returned are the residuals from the linear model: FE ~ log2(sequence region length). See details below for usage recommendations. The default is \code{FALSE}.
  }
  \item{region}{A character vector specifying the sequence region(s) to be analyzed. This can be \code{"UTR5"}, \code{"UTR3"}, \code{"CDS"}, or any combination of the regions.
  }
  \item{comparisons}{A list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the folding energies for the specified sequence region(s) will be returned, and statistical comparisons between gene sets will not be performed. The default is \code{TRUE}. 
  }
  \item{plotType}{If \code{plotOut} is \code{TRUE}, provide a string to indicate the type of output plot to generate. Options are \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"ecdf"}.
  }
  \item{pdfName}{Name to be appended to output PDF files. The default is \code{NULL}.
  }
}

\details{Pre-calculated folding energies supplied with the package are currently available for human and mouse RefSeq releases "rel_109.20201120", and "rel_109.20200923". These values were calculated for all sequence regions using the mfold algorithm, available here: \url{https://www.unafold.org/}.

If the output folding energies will be used in downstream analysis with \link[postNet]{featureIntegration}, consideration should be given to the \code{residFE} argument. If the length of the sequence regions will be included in the \code{featureIntegration} model, it is not recommended to correct the folding energies to the sequence length. However, if length will not be included in \code{featureIntegration}, the correction can be performed.   
}

\value{The output will be a named list of vectors corresponding to each sequence region, with the folding energy for each gene. This list can be used as input for the downstream \link[postNet]{featureIntegration} analysis. The \code{foldingEnergyAnalysis} function can also return PDF files of output plots with statistical comparisons between gene sets of interest, or against background.
}

\references{
If you use the \code{foldingEnergyAnalysis} function in your analysis, with the \code{sourceFE = "load"}, please cite: \cr

M. Zuker. Mfold web server for nucleic acid folding and hybridization prediction. Nucleic Acids Res. 31 (13), 3406-3415, 2003. https://doi.org/10.1093/nar/gkg595
}

\seealso{
\url{https://www.unafold.org/}
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
 

# Compare the folding energy of the 5'UTR for each gene between translationUp genes vs. background, 
# translationDown genes vs. background, and translationUp vs. translationDown genes. 
 FE <- foldingEnergyAnalysis(ptn = ptn,
                               region=c("UTR5"),
                               comparisons = list(c(0,1),c(0,2),c(1,2)),
                               residFE = FALSE, 
                               plotType = "ecdf",
                               sourceFE = "load",
                               plotOut = TRUE,
                               pdfName = "example")
                      
 str(FE)
}
