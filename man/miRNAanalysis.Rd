\name{miRNAanalysis}
\alias{miRNAanalysis}
\title{Perform miRNA target enrichment analysis using a \code{postNetData} object 
}

\description{The \code{miRNAanalysis} function employs \link[gage]{gage} in combination with lists of mRNA targets of miRNA from the 
TargetScan database to identify enrichments in miRNA targets within the gene sets of interest contained in a \code{postNetData} object. 
This function does not identify or predict miRNA binding sites, but rather assesses whether gene sets of interest may be regulated 
by particular miRNAs.     
}

\usage{
miRNAanalysis(ptn,
              miRNATargetScanFile,
              genesSlopeFiltOut = NULL,
              contextScore = -0.2,
              Pct = 0,
              maxSize = 500,
              minSize = 10)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
  \item{miRNATargetScanFile}{The name and path of a TAB delimited text file with miRNA target predictions. The target prediction file must contain
  the headings 'Cumulative weighted context++ score', 'Aggregate PCT', 'Gene Symbol', and 'Representative miRNA', and can be downloaded from the 
  Target Scan database \url{https://www.targetscan.org}. Importantly, the file must be subset to contain predictions for 
  only the desired species prior to running the analysis.
  }
  \item{genesSlopeFiltOut}{If using an \code{Anota2seqDataSet} with analysis results for the "translation" or
  "buffering" regulatory modes, optionally provide a character vector with gene IDs to be excluded from GSEA that 
  is the output of the \link[postNet]{slopeFilt} function. Although optional, this step is strongly 
  recommended. This filtering is not necessary if using "mRNAAbundance", "totalmRNA", or "translatedmRNA". 
  Default is \code{NULL}.
  }
  \item{contextScore}{A negative numeric value specifying the threshold of the cumulative weighted context++ score on which to filter
  miRNA target predictions. Larger negative values of \code{contextScore} indicate stronger biochemical response to a given miRNA (more repression). 
  See details below. Default is \code{-0.2}.
  }
  \item{Pct}{A numeric value between 0 and 1 specifying the threshold of the aggregate Pct score on which to filter
  miRNA target predictions. Larger values of \code{Pct} (closer to 1) indicate a greater probability that sites for a given gene have been
  evolutionarily conserved due to maintenance of miRNA tartgeting (biological function). Pct and cumulative weighted context++ score represent independent
  alternative scoring schemes. Therefore, filtering can be performed on both metrics, or according to one metric at a time.
  See details below. Default is \code{0}.
  }
  \item{maxSize}{Integer specifying the maximal size of a gene set to test. miRNA with more predicted target genes will 
  be exluded. Default is \code{500}.
  }
  \item{minSize}{Integer specifying the minimal size of a gene set to test. miRNAs with fewer predicted target genes will 
  be exluded. Default is \code{10}.
  }
}

\details{For each miRNA-mRNA targeting prediction, two metrics are obtained from the TargetScan database; the cumulative weighted context++ score (CWCS) (see
\url{https://www.targetscan.org/vert_71/docs/context_score_totals.html} for details) and the Pct score (see \url{https://www.targetscan.org/docs/pct.html} for details). 
In order to limit the analysis to high-confidence targeting predictions, filtering can be performed using one, or both of these metrics. Subsequently, \link[gage]{gage}
is applied to identify enrichments in miRNAs predicted to target gene sets of interest based on ranked genes (ranked by the regulatory effect measurement).

Depending on the biological question, it is important to consider which thresholds to use to filter target predictions. Context scores provide a prediction of targeting efficiency, whereas Pct is an estimate of the probability that the miRNA-mRNA targeting interaction is evolutionarily conserved suggesting important biological function. See \url{https://www.targetscan.org/faqs.html} for further details regarding target prediction. 

Note that enrichments in miRNA predicted to target genes that are upregulated (e.g., if the regulatory effect measurement is log2 fold change) that appear
in the results table labelled "greater" can be interpretted as those miRNA that may be dwonregulated or otherwise not active in the experimental condition. Likewise, 
enrichments in miRNA predicted to target genes that are downregulated that appear in the results table labelled "less" can be interpretted as those miRNA that may be
upregulated or active in the experimental condition. 

After running GAGE miRNA target enrichment analysis, the results can be retrieved from the \code{postNetData} object using the S4 method \link[postNet]{ptn_miRNA_analysis} function. To include specific miRNA in downstream \link[postNet]{featureIntegration} analysis, lists of genes targeted by particular miRNA can be retrieved using the S4 method \link[postNet]{ptn_miRNA_to_gene} function. These gene lists can then be included as signatures in the models. 

Note that if the results of an \pkg{anota2seq} analysis are being used for GAGE miRNA target enrichment analysis, it is strongly recommended to first perform filtering using the \link[postNet]{slopeFilt} function.
}

\value{The results returned by the \code{miRNAanalysis} function are compiled into an S4 object of class \code{postNetmiRNA} and stored in the "miRNA" slot of the 
\code{postNetData} object.  The "miRNA_analysis" slot stores a named list with three elements ("greater", "less" and "stats"). As described in the \link[gage]{gage} package, the "stats" element contains the test statistics. The "greater" and "less" elements are data.frames with the results of a two-directional test, each having the 
following columns:

 \itemize{ 
  \item \strong{p.geomean}: The geometric mean of the individual p-values from multiple single array-based gene set 
  tests.
  \item \strong{stat.mean}: The mean of the individual statistics from multple single array-based gene set tests. The
  value denotes the magnitude of the gene-set level changes, and the sign denotes the direction of the changes.
  \item \strong{p.val}: The global p-value or summary of the individual p-values from multiple single array-based gene
  set tests.
  \item \strong{q.val}: The BH adjusted p-value, as implemented by the \pkg{multtest} package.
  \item \strong{set.size}: The number of genes included in the gene set.
  }
The contents of the "miRNA_analysis" slot can be accessed using the \link[postNet]{ptn_miRNA_analysis} S4 method. 
\cr
The "miRNA_to_gene" slot stores a list where each element is a vector of genes that are predicted to be targeted by the miRNA that passed the filtering thresholds for \code{contextScore} and/or \code{Pct}. 

The contents of the "miRNA_to_gene" slot can be accessed using the \link[postNet]{ptn_miRNA_to_gene} S4 method. 
}

\references{

If you use the \code{miRNAanalysis} function in your analysis, please cite:

McGeary SE, Lin KS, Shi CY, Pham T, Bisaria N, Kelley GM, Bartel DP. The biochemical basis of microRNA targeting efficacy. Science Dec 5, (2019). 

Agarwal V, Bell GW, Nam J, Bartel DP. Predicting effective microRNA target sites in mammalian mRNAs. eLife, 4:e05005, (2015). eLife Lens view.

Luo, W., Friedman, M., Shedden K., Hankenson, K. and Woolf, P GAGE: Generally Applicable Gene Set Enrichment for Pathways Analysis. BMC Bioinformatics 2009, 10:161.
}

\seealso{
\link[gage]{gage} \cr
\link[postNet]{ptn_miRNA_to_gene} \cr
\link[postNet]{ptn_miRNA_analysis} \cr
}

\examples{

\dontrun{
# ------------------------------------------------------------------
# Example 1: Running GAGE miRNA target enrichment analysis with gene lists 
# ------------------------------------------------------------------

# load example data:
 data("postNetExample", package = "postNetParcel")
 ptn <- postNetExample$ptn

# Run GAGE miRNA target enrichment analysis:
 ptn <- miRNAanalysis(ptn = ptn,
                     miRNATargetScanFile = 'miRNA_predictions_hsa.txt',
                     contextScore = -0.2,
                     Pct = 0.4)

# Extract the signficant enrichment results from the postNetData object:
 miRNAout <- ptn_miRNA_analysis(ptn = ptn,
                    direction = "less", 
                    threshold = 1)
 str(miRNAout)

# ------------------------------------------------------------------
# Example 2: Running analysis requiring anota2seq analysis slope filtering
# ------------------------------------------------------------------

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
    
# Initialize the postNetData object, where "ads" is an Anota2seqDataSet object
# resulting from an anota2seqRun call:
  ptn <- postNetStart(
      ads = ads,
      regulation = c("translationUp","translationDown"),
      contrast = c(1,1),
      regulationGen = "translation",
      contrastSel = 1,
      selection = "random",
      setSeed = 123,
      source = "load",
      species = "human"
    )
    
# Run slope filtering to remove the genes with unrealistic slopes:
 filtOutGenes <- slopeFilt(ads = ads, 
                          regulationGen="translation", 
                          contrastSel=1)
                          
# Run GAGE miRNA target enrichment analysis:
 ptn <- miRNAanalysis(ptn,
                     genesSlopeFiltOut = filtOutGenes, 
                     miRNATargetScanFile = 'miRNA_predictions_hsa.txt',
                     Pct = 0.7)
  }
}
