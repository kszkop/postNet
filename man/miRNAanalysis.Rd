\name{miRNAanalysis}
\alias{miRNAanalysis}
\title{Perform miRNA target enrichment analysis using a \code{postNetData} object 
}

\description{The \code{miRNAanalysis} function employs \link[gage]{gage} in combination with lists of mRNA targets of miRNA from the TargetScan database to identify enrichments in miRNA targets within the gene sets of interest contained in a \code{postNetData} object. This function does not identify or predict miRNA binding sites, but rather assesses whether gene sets of interest may be regulated by particular miRNAs.    
}

\usage{
miRNAanalysis(ptn,
              miRNATargetScanFile,
              genesSlopeFiltOut = NULL,
              maxSize = 500,
              minSize = 10)
}

\arguments{
  \item{ptn}{a \code{postNetData} object.
  }
  \item{miRNATargetScanFile}{name of a file including a path with predicted (conserved) targets of conserved miRNA families with 
  context++ scores downloaded from the Target Scan database (https://www.targetscan.org). Important, please subset the file for 
  only desired species.
  }
  \item{genesSlopeFiltOut}{If using an \code{Anota2seqDataSet} with analysis results for the "translation" or
  "buffering" regulatory modes, optionally provide a character vector with gene IDs to be excluded from GSEA that 
  is the output of the \link[postNet]{slopeFilt} function. Although optional, this step is strongly 
  recommended. This filtering is not necessary if using "mRNAAbundance", "totalmRNA", or "translatedmRNA". 
  Default is \code{NULL}.
  }
  \item{maxSize}{Integer specifying the maximal size of a gene set to test. miRNAs with more target genes will 
  be exluded. Default is \code{500}.
  }
  \item{minSize}{Integer specifying the minimal size of a gene set to test. miRNAs with fewer target genes will 
  be exluded. Default is \code{10}.
  }
}

\details{For each miRNA-mRNA target pairing, the cumulative weighted context score (CWCS) and the predicted repression magnitude is obtained. Data is filtered out based on that score (the lower the score the better repression) and Site.type  - 1 to 4 with 1 being the best match (i.e. A in position 1 followed by 2-8 with matching nucleotides). Subsequently, gage (https://bioconductor.org/packages/release/bioc/html/gage.html) is applied to identify enriched miRNAs based on ranked genes (ranked by log10 rvm p-value, with direction based on log2 fold change values).
}

\value{
}
\references{

If you use "runStreme()" in your analysis, please cite:
Luo, W., Friedman, M., Shedden K., Hankenson, K. and Woolf, P GAGE: Generally Applicable Gene Set Enrichment for Pathways Analysis. BMC Bioinformatics, 2009, 10:161
}

\seealso{
\link[gage]{gage} \cr
}

\examples{
#First prepare anota2seqUtils class object
# from anota2seq object
a2sU <- anota2seqUtilsStart(ads = anota2seqOutput, regulation =  c("translationUp","translationDown"), contrast = c(1,1), regulationGen = "translation", contrastSel = 1, selection ="random", source = 'load', species = 'human')
# or from geneList
a2sU <- anota2seqUtilsStart(geneList = inList, geneListcolours=c("#F2A104","#00743F"), effectMeasure = quantIn, selection ="random")

# Optionally run slope filtering 
filtOutGenes <- slopeFilt(ads = anota2seqOutput, regulationGen="translation", contrastSel=1, minSlope=-1, maxSlope=2)

# Download file with predicted (conserved) targets of conserved miRNA families with context++ scores from Target Scan database (https://www.targetscan.org) 
# Important, please subset the file for only desired specie:
{downloaded file}_filt <- {downloaded file}[{downloaded file}$Gene.Tax.ID == "9606",]

a2sU <- miRNAanalysis(a2sU = a2sU, miRNATargetScanFile= {downloaded file}_filt , genesSlopeFiltOut = filtOutGenes)

#To extract results, with direction to choose 'greater' or 'less'
resOut <- get_miRNA_results(a2sU = a2sU, direction)

#To extract miRNA targets with miRNAs vector of miRNAs id (for example: c("ptr-miR-497","ptr-miR-16"))
targetsOut <- get_miRNA_geneTargets(a2sU = a2sU, miRNAs)
}
