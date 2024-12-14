# postNet
### Post-transcriptional Network Modeling 
<br>

## Table of contents

* [Overview](https://gitlab.com/krzysztof.szkop/postNet#overview)
* [Citing postNet](https://gitlab.com/krzysztof.szkop/postNet#citing-anota2seqUtils)
* [License](https://gitlab.com/krzysztof.szkop/postNet#license)
* [Before starting](https://gitlab.com/krzysztof.szkop/postNet#before-starting)
	- [Dependencies](https://gitlab.com/krzysztof.szkop/postNet#dependencies)
	- [Installation](https://gitlab.com/krzysztof.szkop/postNet#installation)
	- [Loading](https://gitlab.com/krzysztof.szkop/postNet#loading)
	- [GettingHelp](https://gitlab.com/krzysztof.szkop/postNet#getting-help)
* [Usage](https://gitlab.com/krzysztof.szkop/postNet#usage)
* [Contacts](https://gitlab.com/krzysztof.szkop/postNet#contacts)

------------------------------------------------------------------------

## Overview

------------------------------------------------------------------------

## Citing postNet

Please cite the following article when using __postNet__:

------------------------------------------------------------------------

## License

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg?maxAge=2592000?style=flat)](https://opensource.org/licenses/MIT)

------------------------------------------------------------------------

## Before starting

### Dependencies

postNet requires:

 MEME Suite version >= 5.3.3 (required for motif discovery)

 Cytoscape version 3.8.2 with ClueGo app version >= 2.5.7 (required for network analysis)

 mfold version 3.6 (required for fold energy calculations)

 R version >= 4.1.0 and the following R packages:

* CRAN
    + gridExtra (>= 2.3)
    + dplyr (>= 1.0.7)
    + data.table (>= 1.14.6)
    + curl (>= 4.3.1)
    + seqinr (>= 4.2.8)
    + R.utils (>= 2.10.1)
    + reshape2 (>= 1.4.4)
    + vioplot (>= 0.3.6)
    + stringr (>= 1.4.0)
    + plotrix (>= 3.8.1)
    + gplots (>= 3.1.1)
    + ggplot2 (>= 3.3.4)
    + ggrepel (>= 0.9.1)
    + WriteXLS (> 6.3.0)
    + randomForest (>= 4.7-1.1)
    + igraph (>= 1.3.5)
    + Boruta (>= 7.0.0)
    + ROCR (>= 1.0-11)
    + caret (>= 6.0-93)
    + msigdb (>= 1.6.0)
    + RColorBrewer (>= 1.1-3)

* Bioconductor
    + anota2seq (>= 1.14.0)
    + memes (>= 1.0.0)
    + pqsfinder (>= 2.8.0)
    + GenomicRanges (>= 1.44.0)
    + GOstats (>=2.58.0)
    + Biostrings (>= 2.60.2)
    + ExperimentHub (>= 2.6.0)
    + AnnotationHub (>= 3.6.0)
    + GSEABase (>= 1.60.0)
    + fgsea (>= 1.24.0)
    + org.Hs.eg.db (>= 3.16.0)
    + org.Mm.eg.db (>= 3.16.0)



### Installation

 To install __postNet__ directly from GitLab the *devtools* package is required. If not already installed on your system, run
    
    install.packages("devtools")
	
 Install __postNet__ by
	
	devtools::install_git("https://gitlab.com/krzysztof.szkop/postNet.git",credentials = git2r::cred_user_pass(<User name>, <password>), build_vignettes = TRUE)
    
### Loading

 To load __postNet__ run

	library(postNet)

### Getting help

 Next sections illustrate how to make use of __postNet__ by introducing all functions included in the package and reporting most of the data structures and graphical outputs generated with the default options. Similar information are reported in the vignette returned by
 
	browseVignettes("postNet")
 
 For additional examples and further details about the meaning and usage of all parameters in a function run
 
	?function_name

 or

	help(package = postNet)
 
 A complete reference manual is available [here](https://gitlab.com/krzysztof.szkop/postNet/main/ReferenceManual.pdf).   

 Bugs and errors can be reported at the [issues](https://gitlab.com/krzysztof.szkop/postNet/issues) page on gitlab. Before filing new issues, please read the documentation and take a look at currently open and already closed discussions.

------------------------------------------------------------------------    
    
## Usage - default run

-------------------------------------------------

###Load postNet

    library(postNet)

###Load example of anota2seq (it is one comparison from osmosis)

    data(postNetData)

###Preparation step

####if from anota2seqObject

    ptn <- postNetStart(ads = anota2seqObject, regulation =  c("translationUp","translationDown"), contrast = c(1,1), regulationGen = "translation", contrastSel = 1, region = c("UTR5","CDS","UTR3"), selection ="random", source = 'load', species = 'human')

####Optional slope filtering - it is important if GO/GSEA/GAGE analysis is performed as these genes are excluded from the anota2seq analysis

    filtOutGenes <- slopeFilt(ads = ads, regulationGen="translation", contrastSel = 1, minSlope=-1, maxSlope=2)

####if from gene list

    ptn <- postNetStart(geneList = geneList[1:2], geneListcolours=c("#F2A104","#00743F"), effectMeasure = quantIn, region = c("UTR5","CDS","UTR3"), selection ="random", source = 'load', species = 'human')


###GSEA analysis

####Perform gsea analysis

    ptn <- gseaAnalysis(ptn = ptn, genesSlopeFiltOut = filtOutGenes, collection = 'h')
    
####Extract results

    gseaOut <- ptn_gsea(ptn, threshold = 0.05)
    
####Plot selected terms

    gseaPlot(ptn = ptn,termNames = gseaOut$Term[1:2], genesSlopeFiltOut = filtOutGenes, gseaParam = 1, ticksSize = 0.3)

###GO analysis

####Perform go analysis

    ptn <- goAnalysis(ptn = ptn,genesSlopeFiltOut = filtOutGenes, category=c('BP', 'CC', 'MF', 'KEGG'))
    
####Extract results for the desired category

    BP_GO <- ptn_GO(ptn, category = 'BP', geneList = 'translationUp_c1', threshold = 0.05)
    
####Plot results

    goDotplot(ptn=ptn, category = 'BP', nCategories=10, pool=F)

###GAGE analysis

####Perform go analysis

    ptn <- gageAnalysis(ptn, genesSlopeFiltOut = filtOutGenes, category=c('BP', 'CC', 'MF', 'KEGG'))

####Extract results

    BP_GAGE <- ptn_gage(ptn,category = 'BP',direction = 'greater',threshold = 1)

###miRNA analysis 

####Perform miRNA analysis

    ptn <- miRNAanalysis(ptn,genesSlopeFiltOut = filtOutGenes, miRNATargetScanFile = 'Predicted_Targets_Context_Scores_human.txt')

####Extract results

    miRNAOut <- ptn_miRNA(ptn, direction = 'greater',threshold = 1)          

###Calculate features

####length

    len <- lengthAnalysis(ptn = ptn, region = c("UTR5","CDS","UTR3"), comparisons = list(c(0,1),c(0,2),c(1,2)))
    
####nucleotide content

    content <- contentAnalysis(ptn = ptn, region = c("UTR5","CDS","UTR3"), comparisons = list(c(0,1),c(0,2),c(1,2)), contentIn = c('C','G'))

####uORF presence

    uORFs <- uorfAnalysis(ptn = ptn,comparisons = list(c(0,1),c(0,2),c(1,2)))

####In order to find de novo motifs, run motif analysis

    ptn <- motifAnalysis(ptn = ptn, memePath="/opt/local/bin/", region = c('UTR5'))
    
####motif presence

    motifs <- contentMotifs(ptn = ptn, motifsIn = ptn_motifs(ptn, region = 'UTR5'),region = 'UTR5', comparisons = list(c(1,2)))

####folding energy

    feOut <- foldingEnergyAnalysis(ptn = ptn, region=c("UTR5","CDS","UTR3"), comparisons = list(c(0,1),c(0,2),c(1,2)), residFE = TRUE, plotType = 'ecdf')

####run codonUsage analysis first

    ptn <- codonUsage(ptn = ptn,analysis='codon',comparisons = list(c(1,2)),annotType='ccds',sourceSeq='load',pAdj=0.01,plotHeatmap=TRUE, codonN=1)
    
####calculate codons

    selCodonOut  <- codonCalc(ptn = ptn,analysis='codon', featsel= ptn_codonsSel(ptn, comparison = 1),unit='freq', comparisons = list(c(1,2)))

####calculte signatures

    data(humanSignatures)
    sign <- signCalc(ptn=ptn, signatures = humanSignatures)

###Plot signatures

    signatureFunction(signatureList = humanSignatures[1:2] ,generalName = 'Lee_etal_2014_senescence', dataName= 'Treatment vs Control',colours = c(#FC9272", "#99000D"), ads = anota2seqObject, contrast = 1,pdfName = 'test',xlim = c(-3,3),scatterXY = 5,tableCex=1)

###create input object with calculated features

    features <- c(len[c(1)], content, uORFs, motifs, sign, selCodonOut)

####run feature integration

    featureIntegration(ptn=ptn, features = features, pdfName = 'Test', regOnly=TRUE, allFeat = F, analysis_type ='lm',comparisons = list(c(1,2)))

####Btw there are two new arguments for feature integration function that allow to categorise features into groups . It is still quite experimental and not sure how well it will work in different datasets. Grouping: (vector of respective group for features. I hope you can use names of groups instead of numbers
####lmfeatGroup = c(rep(1,8),rep(2,18),rep(3,13),rep(1,2))
####And colours for each group:
####lmfeatGroupColour = c('#F2C6ED', '#C9E4DE','#FAEDCB’)


------------------------------------------------------------------------

## Contacts

krzysztof.szkop@gmail.com

