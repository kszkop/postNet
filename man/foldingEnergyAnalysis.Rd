\name{foldingEnergyAnalysis}
\alias{foldingEnergyAnalysis}
\title{Compare mRNA folding energy between gene sets  
}

\description{Calculation of folding energy is performed. mfold algorithm is used for calculations (if the function is 
used for calulcuating folding energies (sourceFE = "create"), it requires installed mfold with executables in the path). 
Mfold can be run from provided a fasta file or from prepared by anota2seqUtilsStart function. As the process of calculating 
folding energy is very long (even weeks for a large number of long sequences), it is also possible just to use already 
precalculated values for available reference sequences (available for mouse and human). Furthermore, differences in 
folding energy between selected modes of regulation can be also assessed and plotted.
}

\usage{
foldingEnergyAnalysis(ptn, 
                      sourceFE = "load", 
                      customFileFE = NULL, 
                      residFE = FALSE, 
                      region = NULL, 
                      comparisons = NULL, 
                      plotOut = TRUE, 
                      plotType = "ecdf", 
                      pdfName = NULL)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
  \item{sourceFE}{selecting the source of the folding energy calculations. If "load", it will use precalculated folding 
  energies for selected species ("human" and "mouse" available at the moment). However, there is an option to automatically
  calculate new ("create") based on the selection of regions or fromFasta ( fasta file with sequences and transcript ids as
  names) or load custom 'custom' file with precalculated folding energies calculations in the format: transcript ID, 
  folding energy, length of the sequence. The default is "load". 
  }
  \item{customFileFE}{Name of a file including a path to the customFile. It is required for sourceFE=='custom' or if 
  sourceFE=='create' and fromFasta==TRUE. In the former case, it must be the file with calculated folding energies in 
  the format: transcript ID, folding energy, length of the sequence, or in the latter fasta files with sequences and 
  transcript ids as names. 
  }
  \item{residFE}{Residuals after correction for the length of the sequences should be used instead. The default is TRUE.
  }
  \item{region}{A character vector specifying the sequence region(s) of interest. This can be \code{"UTR5"},
  \code{"UTR3"}, \code{"CDS"}, or any combination of the regions.
  }
  \item{comparisons}{A list of numeric vectors indicating which gene sets in \code{regulation} or 
  \code{geneList} should be used in pair-wise comparisons. The numeric vectors should correspond to the
  indexes of the gene sets, where 0 denotes the background. For example, \code{list(c(0,2),c(1,2))} 
  indicates that the second gene set in \code{regulation} or \code{geneList} will be compared to the
  background (0), and the first and second gene sets will be compared to each other. The default is
  \code{NULL}.
  }
  \item{plotOut}{Logical indicating whether PDF files of plots are generated. If \code{FALSE}, only the 
  fold energies for the sepcified sequence region(s) will be returned, and statistical comparisons between gene 
  sets will not be performed. The default is \code{TRUE}. 
  }
  \item{plotType}{If \code{plotOut} is \code{TRUE}, indicate the type of output plot to generate. Options 
  are, \code{"boxplot"}, \code{"violin"}, or \code{"ecdf"}. The default is \code{"ecdf"}.
  }
  \item{pdfName}{Name to be appended to output PDF files.
  }
}

\details{The function uses mfold algorithm for folding energy calculations (https://www.unafold.org/mfold/software/download-mfold.php).
}

\value{Vector with calculated folding energy values for each gene. The function will also output plots with 
desired comparisons.
}

\references{
# Citation
M. Zuker. Mfold web server for nucleic acid folding and hybridization prediction. 
Nucleic Acids Res. 31 (13), 3406-3415, 2003. https://doi.org/10.1093/nar/gkg595
}

\seealso{
}

\examples{
#First prepare anota2seqUtils class object
# from anota2seq object
a2sU <- anota2seqUtilsStart(ads = anota2seqOutput, regulation =  c("translationUp","translationDown"), contrast = c(1,1), regulationGen = "translation", contrastSel = 1, selection ="random", source = 'load', species = 'human')
# or from geneList
a2sU <- anota2seqUtilsStart(geneList = inList, geneListcolours=c("#F2A104","#00743F"), effectMeasure = quantIn, selection ="random")

feOut <- foldingEnergyAnalysis(a2sU = a2sU, sourceFE = "load", residFE = TRUE, region = c("UTR5","CDS","UTR3"),comparisons = list(c(1,2)),plotOut = TRUE, pdfName = NULL )

#If only calculating folding energies
foldingEnergyAnalysis(a2sU = a2sU, sourceFE = "create", region = c("UTR5"))
}
