\name{featureIntegration}
\alias{featureIntegration}
\title{Feature integration with post-transcriptional network modelling
}

\description{The \code{featureIntegration} function is used to model changes in post-transcriptional regulation and identify cis and/or trans factors (features of mRNA molecules, signatures of upstream regulators, etc.) that can explain the observed regulatory effects.  The method integrates mRNA features and signatures obtained in previous steps of the postNet workflow (or from custom sources), and has several options for implementation, including:
    \itemize{
      \item Stepwise linear regression modelling and network analysis. This analysis uses heirarchical multiple regression to identify features explaining changes in the regulatory effect, rank them according to their importance, and reveals independent and combinatorial effects (for example, co-occurrance of regulatory features in the same mRNA molecules). The results of this analysis can be visualized using network plots showing the association of features to regulatory effects, and the relationships between features.  
     
      \item Random Forest feature selection and classifcation. This analysis implements Random Forest classification using \code{\link[randomForest]{randomForest}}. Feature selection is performed using \code{\link[Boruta]{Boruta}} to identify the set of features that best classify genes according to their regulation. The Random Forest model can then also be used to predict regulation for other gene lists or datasets using the \code{\link[postNet]{rfPred}} function.     
    }
    
Both implementations of \code{featureIntegration} produce statistics and plots visualizating the relationship between the regulatory effect and the selected features. The relationships between selected features identified by \code{featureIntegration} and regulatory effects can be further explored using UMAP visualizations with the \code{\link[postNet]{plotFeaturesMap}} function.    
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
  \item{regOnly}{If \code{analysis_type = "lm"}, provide a logical specifying whether modelling should be performed using only regulated genes (i.e., those included in \code{geneList}), or with all genes in the dataset. If \code{regOnly = TRUE}, only regulated genes will be used, if \code{regOnly = FALSE}, all genes will be considered. Default is \code{TRUE}. 
}
  \item{allFeat}{If \code{analysis_type = "lm"}, provide a logical specifying which input features should be including in the stepwise linear regression modelling. If \code{allFeat = TRUE}, all input features will be included, if \code{FALSE}, only features passing the FDR threshold specified with \code{fdrUni} in univariate models will be included. Default is \code{FALSE}.
}
  \item{useCorel}{If \code{analysis_type = "lm"}, provide a logical specifying if substantial correlations between features should be displayed as edges in network plots and coefficients be returned along with the final model. If \code{TRUE}, Pearson's correlation coefficients will be calculated between associated features, allowing visualization of how features co-vary in network plots. A threshold for the strength of association between features in order for correlation coefficients to be calculated must be provided using the \code{covarFilt} parameter. If \code{FALSE}, edges in metwork plots will instead represent the magnitude of the change in F-value for a given feature between sequential steps of stepwise regression to display the strength of associations between features, and no correlation coefficients will be calculated. Default is \code{TRUE}.    
}
  \item{covarFilt}{If \code{analysis_type = "lm"}, provide an integer specifying the minimum threshold for the change in F-value between sequential steps of stepwise regression to display the strength of associations between features in network graphs. Larger changes in F-values upon addition of a new feature to the model are indicative of a stronger association between features. Default is \code{20}, indicating that edges between feature nodes in network plots will be displayed between features where there was a change of 20 or more in the F-value for one feature upon addition of the other to the model. Note that this parameter is only applied to network analysis and not F-table visualizations.   
}
  \item{NetModelSel}{If \code{analysis_type = "lm"}, provide a string specifying how the variance explained should be displayed in the network plot. If \code{NetModelSel = "omnibus"}, the total variance explained by the omnibus model (the smallest model to explain the greated proportion of variance) will be displayed. The size of each node (representing features) is proportional to the variance explained, where covariance is assigned to the most influential feature (a rank-based greedy model). If \code{NetModelSel = "adjusted"}, the covariance will be eliminated and the values displayed will reflect the adjusted total variance explained, and the independent contribution of each feature. Default is \code{"omnibus"}. 
}
  \item{comparisons}{If \code{regOnly = TRUE}, provide a list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}. If both  \code{regulation} and \code{geneList} are provided, \code{geneList} is appended to the end so numbering will continue sequentially.
}
  \item{fdrUni}{If \code{analysis_type = "lm"}, provide a numeric value between 0 and 1 specifying the FDR threshold to select significant features in univariate analyses. When \code{allFeat = FALSE}, only the features passing this threshold will be included in the stepwise regression modelling. Default is \code{0.05}.  
}
  \item{stepP}{If \code{analysis_type = "lm"}, provide a numeric value between 0 and 1 specifying the p-value threshold to determine if features should be retained in the omnibus model after each step of stepwise regression. Default is \code{0.05}.
}
  \item{pdfName}{Name to be appended to output PDF files.
}
}

\details{The \code{featureIntegration} analysis assesses whether a catalogue of cis and/or trans regulatory features appear to modulate the observed post-transcriptional regulation. Details of the implementations and considerations for feature inputs are described below.

\strong{Stepwise linear regression modelling and network analysis:} \code{featureIntegration} uses stepwise linear regressions to model changes in post-transcriptional regulation across subsets of mRNAs and assesses the contribution of each feature included in the modelling in a hierarchical manner. This allows features to be ranked according to their ability to explain changes in the observed regulatory effect. Due to the hierarchical approach, identification of both distinct (independent), and overlapping (covarying) associations between features and regulatory effects is possible. \code{featureIntegration} with stepwise linear regression is performed in three phases:

 \itemize{ 
  \item \strong{Phase 1}: First, each candidate feature is evaluated separately in univariate linear models to identify significant associations between the regulatory effect and individual features.
  
  \item \strong{Phase 2}: Second, starting with the feature that best explained changes in the regulatory effect from univariate models, forward stepwise regression is performed by adding features to the model in an iterative fashion, keeping covariance assigned to the most influential feature (a rank-based greedy model). In each step, the best performing model (the feature with the strongest association with the regulatory effect) is retained and features that fail to explain additional variance are discarded. The resulting omnibus model represents the smallest set of features that can explain the greatest proportion of variance in the regulatory effect. This step ranks features according to their importance, and reveals dependence between them. For example, regulation may depend partially on the length of the 5'UTR, but also on the GC content, or the presence of a motif, etc.  
  
  \item \strong{Phase 3}: Finally, the independent contribution of each feature identified in phase 2 is determined. This is done by removing covariance from the omnibus model (that was previously assigned to the most influential features in phase 2) to provide the adjusted contribution of each feature. This phase provides the percentage of variance in the regulatory effect explained by a given feature, independent of all other features.  
  }

The results of all three phases of the stepwise regression modelling are summarized in table outputs and visualized as network plots (described below).

\strong{Random Forest feature selection and classifcation:} \code{featureIntegration} can also apply Breiman's random forest algorithm (implemented in the \code{\link[randomForest]{randomForest}} package) to classify genes according to their regulation, and identify features associated with regulatory classes (for example, translationally activated or suppressed genes). This implementation is carried out in two steps: 

 \itemize{ 
  \item \strong{Pre-modelling and feature selection}: Genes are first divided into training (70\%) and validation (30\%) sets. Modelling is performed and importance is calculated for all features included in the input. Feature selection is then performed using \code{\link[Boruta]{Boruta}} to identify all input features relevant to model performance. 
  
  \item \strong{Final modelling}: Next, modelling using \code{\link[randomForest]{randomForest}} is repeated using the set of selected features. The performance of the model is assessed using Receiver Operating Characteristic (ROC) curves (implemented with \code{\link[ROCR]{ROCR}}), and feature importance from the final model is reported as Mean Decreased Accuracy (see \code{\link[randomForest]{importance}} for details), where higher values indicate more important features.   
  }

\strong{Feature inputs for modelling:} Many different types of features can be included in modelling with \code{featureIntegration}. From the postNet workflow, the outputs of the \code{\link[postNet]{lengthAnalysis}}, \code{\link[postNet]{contentAnalysis}}, \code{\link[postNet]{uorfAnalysis}}, \code{\link[postNet]{foldingEnergyAnalysis}}, \code{\link[postNet]{contentMotifs}}, and \code{\link[postNet]{codonCalc}} functions can be supplied directly as features. Gene signatures, such as those included with the package (see \code{\link[postNet]{humanSignatures}} or \code{\link[postNet]{mouseSignatures}}) can also be converted to feature inputs using the \code{\link[postNet]{signCalc}} function. In addition, custom feature inputs can be supplied depending on the application and research question.

In general, features that can be included in modelling can be either numeric, or categorical. Numeric features can be both continuous (for example, sequence region nucleotide content or length), or discrete (for example, the number of uORFs in 5'UTRs). Categorical features can also be included by converting to binary variables. This is useful for evaluating gene signatures in relation to regulatory outcomes, but can also be applied in many scenarios where mRNA can be divided into two classes (for example those with, or without a certain property, etc.). It is also possible to supply categorical features with a directionality. For example, a feature that can be increased, decreased, or unchanged could be coded as 1, -1, or 0 for each gene to be included in modelling.    

Careful considation should be given when selecting the catalogue of features to include in modelling. For example, highly, or perfectly correlated variables cannot be supplied together in the same model. This can be relevant when examining the percentage of nucleotide content for each base, since the content of G and C may often be highly correlated, and the sum of A, T, G, and C will add to 100\%. In instances where features are too highly correlated, and warning message will be displayed and one of the correlated features must be removed from the input.  In addition, when enumerating some features such as codon composition or folding energy, these values can be corrected for the length of the sequence region. However, it would then not be advisable to also include the sequence region length as a separate feature in the modelling. Please see the postNet vignette for more in-depth examples and discussion on the selection of features and interpretation of results.  
}

\value{The \code{featureIntegration} function returns an updated \code{postNetData} object with results compiled and stored in the \code{features} and \code{featureIntegration} slots.

The \code{features} slot stores a \code{data.frame} with the features used as input for feature integration modelling, where rows correspond to genes, and columns correspond to features. The input feature \code{data.frame} can be retrieved using the \link[postNet]{ptn_features} S4 method.   

The stepwise regression implementation of \code{featureIntegration} will store results for each comparison in the \code{lm} slot. This includes the results of univariate and stepwise regression modelling, and the final omnibus and adjusted models. Modelling results can be retrieved using the \link[postNet]{ptn_model} S4 method. The features selected as significant in the omnibus model are listed in the \code{selectedFeatures} slot, along with the proportion of variance in the regulatory effect explained by each feature (either in omnibus, or adjusted models depending on the selection in the \code{NetModelSel} parameter). Selected features can be retrieved using the \link[postNet]{ptn_selectedFeatures} S4 method. Finally, the network analysis results are stored as an \link[igraph]{igraph} object, and can be retrieved using the \link[postNet]{ptn_networkGraph} S4 method. 

The resulting network plot is output as a PDF. Modelling is summarized in a PDF file with a table displaying results for features significant in univariate models, a table of F-values from linear models in stepwise regression, and a table summarizing features in the final omnibus and adjusted models. Optionally, if \code{useCorel = TRUE}, a table of Pearson's correlation coefficients between features will also be displayed. For F-value tables, the feature explaining the greatest proportion of variance in the regulatory effect at each step is highlighted in green, while those removed from the model due to inability to singificantly explain variance are marked in red. Features highlighted in orange indicate those where there was a substantial change in F-value (an increase or decrease of at least 50\%) following the addition of another feature to the model. Larger changes in F-value are indicative of a stronger relationship between features. Note that this highlighting is independent of the threshold applied with the \code{covarFilt} parameter, which controls edge selection for network analyses. For more details please refer to the postNet vignette.    

The random forest implementation of \code{featureIntegration} will store results for each comparison in the \code{rf} slot. The pre-model, feature selection, and final model are stored as objects of class \code{\link[randomForest]{randomForest}} and \code{\link[Boruta]{Boruta}}. The selected features selected are listed in the \code{selectedFeatures} slot, along with their importance (Mean Decrease Accuracy) from the final model. Feature importance calculated by \code{\link[Boruta]{Boruta}} is plotted is output as a PDF file. PDF files are also output with ROC curves for the final model, and the importance of selected features.  

Both implementations of \code{featureIntegration} output PDF files with sactterplots for each individual feature identified in final stepwise regression or random forest models, with Pearson's correlation coefficient between the feature and regulatory effect (two-sided test). The linear trend line, and cubic smoothing spline (allowing assessment of linearity) are displayed.  
}

\references{
If you use random forest analysis, please cite:

Sing T, Sander O, Beerenwinkel N, Lengauer T (2005). “ROCR: visualizing classifier performance in R.” Bioinformatics, 21(20), 7881. http://rocr.bioinf.mpi-sb.mpg.de.

Liaw A, Wiener M (2002). “Classification and Regression by randomForest.” R News, 2(3), 18-22. https://CRAN.R-project.org/doc/Rnews/. 

Witold R. Rudnicki M (2010). “Feature Selection with the Boruta Package.” Journal of Statistical Software, 36(11), 1–13. https://doi.org/10.18637/jss.v036.i11.
}

\seealso{
\code{\link[Boruta]{Boruta}} \cr
\code{\link[randomForest]{randomForest}} \cr
\code{\link[ROCR]{ROCR}} \cr
\code{\link[postNet]{rfPred}} \cr
\code{\link[postNet]{plotFeaturesMap}} \cr
\code{\link[postNet]{ptn_model}} \cr
\code{\link[postNet]{ptn_features}} \cr
\code{\link[postNet]{ptn_selectedFeatures}} \cr
\code{\link[postNet]{ptn_networkGraph}} \cr
}

\examples{

# load and create example data:
 data(postNetExample)

 myGenes <- postNetExample$geneList
 myBg <- postNetExample$background
 myEffect <- postNetExample$effect

# Initialize the postNetData object with custom gene lists:
 ptn <- postNetStart(
   geneList = myGenes,
   geneListcolours = c("#7FB7BE","#DB7F67"),
   customBg = myBg,
   effectMeasure = myEffect,
   selection = "random",
   source = "load", 
   species = "human"
 )
 
# Prepare the list of pre-calculated features to be used in feature integration modelling 
 myFeatures <- postNetExample$features
 str(myFeatures)

# Run feature integration modelling using stepwise regression
 ptn <- featureIntegration(ptn = ptn, 
                          features = myFeatures, 
                          pdfName = "Example", 
                          regOnly = TRUE, 
                          allFeat = FALSE,
                          analysis_type = "lm",
                          covarFilt = 20,
                          comparisons = list(c(1,2)),
                          NetModelSel = "omnibus")
                          
# Run feature integration modelling using random forest classification
 ptn <- featureIntegration(ptn = ptn, 
                          features = myFeatures, 
                          pdfName = "Example", 
                          analysis_type ="rf",
                          comparisons = list(c(1,2)))
                          
}
