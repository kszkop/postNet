\name{featureIntegration}
\alias{featureIntegration}
\title{Feature integration with post-transcriptional network modelling
}

\description{To enable efficient assessment of features underlying alterations in translational efficiencies we developed postNet, which includes methods to quantify and identify RNA features (e.g. nucleotide content, length of UTRs and CDS, folding energy of UTRs and CDS, RNA motifs, uORFs, codon composition etc.) associating with selective alterations in mRNA translation and statistical analysis of whether the catalogue of identified features appear to modulate translation in a distinct or overlapping fashion.

PostNet uses stepwise linear regressions to model changes in post-transcriptional regulation (translation, or translational offsetting) across subsets of regulated mRNAs and assesses the contribution of each feature included in the modelling in a hierarchical manner. This allows postNet to rank mRNA features according to their ability to explain changes in gene expression regulation. Furthermore, because of the hierarchical approach, the interdependence between different features is revealed. 

In the first phase of postNet modelling, each candidate feature is evaluated separately (i.e. using univariate linear models). This allows assessment of associations between changes in post-transcriptional regulation and each selected feature to provide an overview. Moreover, the univariate output is used to identify of a set of features for further modeling in phase two of the postNet analysis. In the second phase, postNet selects the feature that best explained the regulatory changes from phase one and adds additional features to the model in an iterative fashion, keeping covariance between features in the model assigned to the feature which appeared most influential. During each iteration, the best performing model (i.e. the feature showing the strongest association with the regulatory effect) is kept and features that do not explain any additional variation are discarded. This generates a hierarchical output that reveals dependence between RNA features (i.e. co-occurrence in mRNAs) identified under phase one. In the third and final phase of postNet analysis, the independent contribution of each of the features that emerged as significant from phase two is determined (i.e. in contrast to assigning covariance to the most influential feature under phase 2, such covariance is eliminated from the output to derive the independent contribution of each feature during phase three).
}

\usage{
featureIntegration(ptn, 
                  features, 
                  lmfeatGroup = NULL, 
                  lmfeatGroupColour = NULL, 
                  analysis_type, 
                  regOnly = TRUE, 
                  allFeat = FALSE, 
                  useCorel = TRUE, 
                  covarFilt = 20, 
                  NetModelSel = "omnibus", 
                  comparisons = NULL, 
                  fdrUni = 0.05, 
                  stepP = 0.05, 
                  pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{features}{A named list of numeric vectors corresponding to features that will be included in modelling of regulatory effects. Each vector element in the list must have names corresponding to gene IDs. These vectors quantifying features can be the outputs of previous steps of the analysis, including from \code{contentAnalysis}, \code{lengthAnalysis}, etc.. Custom inputs can also be included if they are named numeric vectors with geneIDs matching those in the \code{postNetData} object. Finally, lists of geneIDs without numeric values, for example genes regulated downstream of a pathway, can be converted to binary feature vectors using the \link[postNet]{signCalc} function.  
}
  \item{lmfeatGroup}{If \code{analysis_type = "lm"}, optionally provide an character vector specifying how \code{features} will be grouped in the network plot. Groups of features will be circled, with circle colours corresponding to the groups specified. For example, if the input list for \code{features} includes 5'UTR length, uORFs, enriched codons, 3'UTR length, and a 3'UTR motif, a grouping vector could be: \code{c("UTR5", "UTR5", "CDS", "UTR3", "UTR3")}, to organize the layout of features in the network plot according to regions of the mRNA sequence. Grouping vectors must be the same length as the input \code{features} list. Default is \code{NULL}. 
}
  \item{lmfeatGroupColour}{If \code{analysis_type = "lm"}, optionally provide a named character vector specifying colours for the different groups defined with \code{lmfeatGroup} that will be used in the network plot. Colour vectors must be the same length as the number of groups in \code{lmfeatGroup}, and be named according to group. For example, \code{c("red", "blue", "green")} with names \code{c("UTR5", "CDS", "UTR3")}. Default is \code{NULL}.
}
  \item{analysis_type}{A string specifying the method for the analysis. The options are \code{"lm"} to perform stepwise linear regression modelling, or \code{"rf"} to perform a random forest-based analysis. 
}
  \item{regOnly}{Logical specifying whether modelling should be performed using only regulated genes (i.e., those included in \code{geneList}), or with all genes in the datase. If \code{regOnly = TRUE}, only regulated genes will be used, if \code{regOnly = FALSE}, all genes will be considered. Default is \code{TRUE}. 
}
  \item{allFeat}{If \code{analysis_type = "lm"}, provide a logical specifying which input features should be including in the stepwise linear regression modelling. If \code{allFeat = TRUE}, all input features will be included, if \code{FALSE}, only features passing the FDR threshold specified with \code{fdrUni} in univariate models will be included. Default is \code{FALSE}.
}
  \item{useCorel}{If \code{analysis_type = "lm"}, provide a logical specifying if substantial correlations between features should be displayed as edges in network plots and coefficients be returned along with the final model. If \code{TRUE}, Pearson's correlation coefficients will be calculated between associated features, allowing visualization of how features co-vary in network plots. A threshold for the strength of association between features in order for correlation coefficients to be calculated must be provided using the \code{covarFilt} parameter. If \code{FALSE}, edges in metwork plots will instead represent the magnitude of the change in F-value for a given feature between sequential steps of stepwise regression to display the strength of associations between features, and no correlation coefficients will be calculated. Default is \code{TRUE}.    
}
  \item{covarFilt}{If \code{analysis_type = "lm"}, provide an integer specifying the minimum threshold for the change in F-value between sequential steps of stepwise regression to display the strength of associations between features in network graphs. Larger changes in F-values upon addition of a new feature to the model are indicative of a stronger association between features. Default is 20, indicating that edges between feature nodes in network plots will be displayed between features where there was a change of 20 or more in the F-value for one feature upon addition of the other to the model.   
}
  \item{NetModelSel}{Values for the variance explained for the network plot. The default is "Omnibus". Aternatively "Adjusted" or "Univariate" can be provided. 
}
  \item{comparisons}{If \code{regOnly = TRUE}, provide a list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}. If both  \code{regulation} and \code{geneList} are provided, \code{geneList} is appended to the end so numbering will continue sequentially.
}
  \item{fdrUni}{
}
  \item{stepP}{
}
  \item{pdfName}{Name to be appended to output PDF files.
}
}

\details{
}

\value{
}

\references{
}

\seealso{
}

\examples{
}
