# anota2seqUtils
### Tools for post-anota2seq analysis 
<br>

## Table of contents

* [Overview](https://gitlab.com/krzysztof.szkop/anota2seqUtils#overview)
* [Citing anota2seqUtils](https://gitlab.com/krzysztof.szkop/anota2seqUtils#citing-anota2seqUtils)
* [License](https://gitlab.com/krzysztof.szkop/anota2seqUtils#license)
* [Before starting](https://gitlab.com/krzysztof.szkop/anota2seqUtils#before-starting)
	- [Dependencies](https://gitlab.com/krzysztof.szkop/anota2seqUtils#dependencies)
	- [Installation](https://gitlab.com/krzysztof.szkop/anota2seqUtils#installation)
	- [Loading](https://gitlab.com/krzysztof.szkop/anota2seqUtils#loading)
	- [GettingHelp](https://gitlab.com/krzysztof.szkop/anota2seqUtils#getting-help)
* [Usage](https://gitlab.com/krzysztof.szkop/anota2seqUtils#usage)
* [Contacts](https://gitlab.com/krzysztof.szkop/anota2seqUtils#contacts)

------------------------------------------------------------------------

## Overview

------------------------------------------------------------------------

## Citing anota2seqUtils

Please cite the following article when using __anota2seqUtils__:

------------------------------------------------------------------------

## License

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg?maxAge=2592000?style=flat)](https://opensource.org/licenses/MIT)

------------------------------------------------------------------------

## Before starting

### Dependencies

anota2seqUtils requires:

 MEME Suite version >= 5.3.3 (required for motif discovery)

 Cytoscape version 3.8.2 with ClueGo app version >= 2.5.7 (required for network analysis)

 mfold version 3.6 (required for fold energy calculations)

 R version >= 4.1.0 and the following R packages:

* CRAN
    + gridExtra (>= 2.3)
    + dplyr (>= 1.0.7)
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

* Bioconductor
    + anota2seq (>= 1.14.0)
    + memes (>= 1.0.0)
    + pqsfinder (>= 2.8.0)
    + GenomicRanges (>= 1.44.0)
    + GOstats (>=2.58.0)

### Installation

 To install __anota2seqUtils__ directly from GitLab the *devtools* package is required. If not already installed on your system, run
    
    install.packages("devtools")
	
 Install __anota2seqUtils__ by
	
	devtools::install_git("https://gitlab.com/krzysztof.szkop/anota2seqUtils.git",credentials = git2r::cred_user_pass(<User name>, <password>), build_vignettes = TRUE)
    
### Loading

 To load __anota2seqUtils__ run

	library(anota2seqUtils)

### Getting help

 Next sections illustrate how to make use of __anota2seqUtils__ by introducing all functions included in the package and reporting most of the data structures and graphical outputs generated with the default options. Similar information are reported in the vignette returned by
 
	browseVignettes("anota2seqUtils")
 
 For additional examples and further details about the meaning and usage of all parameters in a function run
 
	?function_name

 or

	help(package = anota2seqUtils)
 
 A complete reference manual is available [here](https://gitlab.com/krzysztof.szkop/anota2seqUtils/main/ReferenceManual.pdf).   

 Bugs and errors can be reported at the [issues](https://gitlab.com/krzysztof.szkop/anota2seqUtils/issues) page on GitHub. Before filing new issues, please read the documentation and take a look at currently open and already closed discussions.

------------------------------------------------------------------------    
    
## Usage

-------------------------------------------------

####Not sure yet whether it will be separete step 

The function extracts results from anota2seq object:

    anota2seqGetDirectedRegulations(ads)

__Output:__ List of regulated genes for each contrast. Including the background.

###### ads

anota2seq object

-------------------------------------------------

### 1. Functional analysis:


###### 1.1 ClueGo network analysis

This function inputs gene lists from anota2seq objects and/or datasets and performs Gene Ontology analysis using ClueGO plug-in on Cytoscape.

**Please make sure that Cytoscape is launched before using this function.**

    anota2seqtoClueGO(ads,contrast,externalGeneList=NULL,organismSelection,geneList2=NULL,outputFolderName,customReferenceListExt=NULL,regType,Evidences="All",nodeShape="Ellipse",nodeShape2="Ellipse",clusterColor = "#FF5733",clusterColor2 = "#900C3F",minGene = 3,minPercentageMapped = 4,minGene2 =3, minPercentageMapped2= 4, restrictions=FALSE,visualStyle="ShowGroupDifference", selectOntology=NULL,enrichmentType="Enrichment (Right-sided hypergeometric test)",treeLevelsSpecificity=c(3,8),allTreeLevels=FALSE,pCutoff=0.05,termFusion= NULL, multipleTestingCorrection=TRUE, correctionType="Benjamini-Hochberg",midPval=FALSE,useDoubling=FALSE,setKappa=0.4,grouping=TRUE,groupingType="Kappa Score",groupColoring="Highest Significance",initGroupSize=1, PercGroupsMerge=50, PercTermsMerge=50 )

__Output:__ A folder in the working directory that contains the log file and analysis results from ClueGO.

###### ads
anota2seq Object

###### contrast
contrast in anota2seq object in which gene vectors will be analysed (here select only one contrast)

###### externalGeneList
Optional. A vector that contains the gene names.Do NOT use this when you are using Anota2seq object.

###### organismSelection
Organism, for now only "Homo Sapiens" and "Mus Musculus" are Present

###### geneList2
Optional.A second vector that contains gene names, used in order to compare two gene lists. e.g. Compare gene list from anota2seq object or externalGeneList vs gene list from geneList2.

###### outputFolderName
The folder where the results are going to be saved in.

###### customReferenceListExt
Optional. If externalGeneList is inputted, a vector that contains all the gene names that are tested need to be inputted to customReferenceListExt parameter.

###### regType 
Regulation Type such as "translationUp","translationDown","bufferingmRNAUp","bufferingmRNADown","mRNAAbundanceUp","mRNAAbundanceDown". Do not add this in case an external dataset is used.

###### Evidences
Default is "All". The evidences that can be used to show the relationship between the functional term and gene/gene sets. Can be set up as "All","All_Experimental_(EXP,IDA,IPI,IMI,IGI,IEP)","All_without_IEA","EXP (Inferred from Experiment)","IBA (Inferred from Biological Aspect of Ancestor)","IBD (Inferred from Biological Aspect of Descendent)" "IC (Inferred by Curator)", "IDA (Inferred from Direct Assay)", "IEA (Inferred from Electronic Annotation)","IEP (Inferred from Expression Pattern)", "IGC (Inferred from Genomic Context)","IGI (Inferred from Genetic Interaction)", "IKR (Inferred from Key Residues)" , "IMP (Inferred from Mutant Phenotype)", "IPI (Inferred from Physical Interaction)","ISA (Inferred from Sequence Alignment)", "ISM (Inferred from Sequence Model)", "ISO (Inferred from Sequence Ontology)", "ISM (Inferred from Sequence Model)","ISS (Inferred from Sequence of Structural Similarity)", "NAS (Non-traceable Author Statement)","ND (No biological Data available)","NR (Not Recorded)", "RCA (Inferred from Reviewed Computational Analysis)", "TAS (Traceable Author Statement)".

###### nodeShape
Default is "Ellipse".The shape of the nodes on the graph that will be plotted. Can be set as "Ellipse","Diamond","Hexagon","Octagon","Parallelogram","Rectangle","Round Rectangle","Triangle","V". 

###### nodeShape2
Used when a second gene list for comparison is available. Default is "Ellipse".Can be set as "Ellipse","Diamond","Hexagon","Octagon","Parallelogram","Rectangle","Round Rectangle","Triangle","V".

###### clusterColor
Color of groups when comparing two gene sets. Terms associated with the first gene set will,by default, have "#FF5733" as color. 

###### clusterColor2
Color of groups when comparing two gene sets. Terms associated with the second gene set will,by default, have "#900C3F" as color.

###### minGene
Default = 3. Minimum number of genes needed to define a functional term.

###### minPercentageMapped
Default = 4. Minimum percentage of overlap between the inputted gene list and all the genes associated with a term.

###### minGene2
Used when a second gene list is available.See minGene for definition.

###### minPercentageMapped2
Used when a second gene list is available.See minPercentageMapped for definition.

###### restrictions
Default = FALSE. TRUE for no restrictions in number and percentage per term.

###### visualStyle
Default = ShowGroupDifference". Can be "ShowGroupDifference" to color the terms based on groups or "ShowSignificanceDifference" to color the terms based on their significance levels or "ShowClusterDifference" to color the terms based on the gene list that they are coming from. In case this is set ti "ShowClusterDifference" the terms will get the colors defined in ClusterColor parameter.

###### selectOntology
Default = c("3;Ellipse","8;Hexagon") for Homo Sapiens and c("2;Ellipse","7;Hexagon") for Mus musculus. Order of Ontology;Shape for that Ontology."Please select the ontologies and define the shapes.

- For organism "Homo Sapiens":
    * 0.Human-diseases (CAN ONLY BE USED FOR HOMO SAPIENS)
    * 1.CORUM-3.0-FunCat-MIPS 
    * 2.Chromosomal-Locatio 
    * 3.BiologicalProcess-EBI-UniProt-GOA-ACAP-ARA 
    * 4.CellularComponent-EBI-UniProt-GOA-ACAP-ARA 
    * 5.ImmuneSystemProcess-EBI-UniProt-GOA-ACAP-ARA 
    * 6.MolecularFunction-EBI-UniProt-GOA-ACAP-ARA 
    * 7.ProteinDomain 
    * 8.KEG  
    * 9.REACTOME Pathway 
    * 10.REACTOME Reactions
    * 11.WikiPathway 

- For organism "Mus musculus":
    * 0.CORUM-3.0-FunCat-MIP 
    * 1.Chromosomal-Locatio 
    * 2.BiologicalProcess-EBI-UniProt-GOA-ACAP-ARA 
    * 3.CellularComponent-EBI-UniProt-GOA-ACAP-ARA 
    * 4.ImmuneSystemProcess-EBI-UniProt-GOA-ACAP-ARA 
    * 5.MolecularFunction-EBI-UniProt-GOA-ACAP-ARA 
    * 6.ProteinDomain 
    * 7.KEGG
    * 8.REACTOME Pathway 
    * 9.REACTOME Reaction 
    * 10.WikiPathway
 
Available shapes are Ellipse,Diamond,Hexagon,Octagon,Parallelogram,Rectangle,Round Rectangle,Triangle,V The format needs to be as e.g.:'2;Ellipse  If you use more than one Ontologies, make this a character vector."

###### enrichmentType
Default = "Enrichment (Right-sided hypergeometric test)". Can be set to "Enrichment (Right-sided hypergeometric test)" to find enriched GO terms for a given gene set, "Depletion (Left-sided hypergeometric test)" to find depleted GO terms for a given gene set, "Enrichment/Depletion (Two-sided hypergeometric test)" to find enriched and depleted GO terms for a given gene set.

###### treeLevelsSpecificity
Default is c(3,8), A vector of min and max tree levels. Vector composed of two elements.Defines the specificity of the terms. Lower the value lower the specificity.

###### allTreeLevels
Default = FALSE. Setting min max tree levels to define the network sensitivity. If TRUE, it searches on all tree levels. When false, treeLevels are c(3,8).

###### pCutoff
Default = 0.05. The p-value cutoff.

###### TermFusion
Default = NULL.Fusion of GO parent-child terms based on similar associated genes.

###### multipleTestingCorrection
Default= TRUE. Logical that indicates if multiple testing correction is going to be applied.

###### correctionType
Default = "BH". If Multiple Testing Correction set as TRUE, this parameter indicates which test is going to be used for multiple test correction. Possible tests are Benjamini-Hochberg","Bonferroni","Bonferroni step down","None"

###### midPval 
Default= FALSE. When set to TRUE,it performs a less conservative hypergeometric test.

###### useDoubling
Default= FALSE. When set to TRUE,it performs a less conservative two-tailed hypergeometric test.

###### setKappa
Default = 0.4. Used to define term-term interrelations and functional groups based on shared genes between terms,associates the terms into functional groups. Higher the kappa,better the clustering. 

###### grouping
Default = TRUE. Logical that indicates if the terms will be grouped. 

###### groupingType
Default = "Kappa Score". Can be set as "Kappa Score" or "Tree". If Kappa Score is selected, the terms will be grouped based on kappa score whereas if Tree is selected, the terms will be grouped as hierarchical trees. 

###### groupColoring 
Default = "Random", Could be set as "Random" which would color the groups in a random manner or "Fix" to color the groups in a fixed manner.

###### groupLeadingTerm
Default= "Highest Significance". Used to set the criteria that defines a leading term, can be one of the terms as follows: ("Highest Significance","#Genes / Term","%Genes / Term","%Genes / Term vs Cluster")

###### initGroupSize
Default=1

###### PercGroupsMerge
Default=50 Overlap between two gene sets associated with a term to merge the functional groups.

###### PercTermsMerge
Default=50 Overlap between two sets of terms associated a functional to merge the functional groups.

**Please see :** *Bindea G, Mlecnik B, Hackl H, et al. ClueGO: a Cytoscape plug-in to decipher functionally grouped gene ontology and pathway annotation networks. Bioinformatics. 2009;25(8):1091-1093. doi:10.1093/bioinformatics/btp101* and *https://apps.cytoscape.org/apps/cluego* for more detailed information about the parameters.


###### 1.2 GO analysis (GOstat)

###### 1.3 Gage

###### 1.4 Signatures

This function plots each set of genes in each signature on scatter plot and ecdf for each of the regulatory mode.

    signatureFunction(signatureList, ads,  contrast, generalName, dataName, effects_names=c('total mRNA log2FC', 'polysome associated mRNA log2FC','buffering log2FC','translation log2FC'), colours, xlim=NULL, scatterXY=NULL, tableCex, pdfName=NULL)

__Output:__ Plot

###### ads
anota2seq object

###### contrast
contrast in anota2seq object in which gene vectors will be analysed (here select only one contrast)

###### dataName
name of the the comparison used as dataset

###### signatureList
Named list of vectors of genes. Names will be used in plots

###### generalName
general name of the signatures to be used for plotting

###### effects_names
Default: c('total mRNA log2FC', 'polysome associated mRNA log2FC','buffering log2FC','translation log2FC'), used for naming axis on the plots

###### colours
Provide colours for each respective signature

###### xlim
Default: NULL and the scale for axis of ecdfs is automatically calculated. Option provide as standard xlim function for example c(-3,3)

###### scatterXY
Default: NULL and the scale for axis of scatter plot is automatically calculated. Option provide one number representing scale for axis.

###### tableCex
size of the table with quantitative results (similar use as cex function)

###### pdfName
Optional, name of the output pdf file with plot

Alternatively there is an option to plot all the signatures using heatmap.


    signaturesHeatmap(signatureList, ads,  contrast, contrastNames, unit, RegMode)

__Output:__ Plot

###### ads
anota2seq object

###### contrast
contrasts in anota2seq object in which gene vectors will be analysed. It is possible to compare contrats to each other by providing multiple contrasts

###### contrastNames
names of each respective contrasts to be used for plotting

###### signatureList
Named list of vectors of genes. Names will be used for plotting

###### unit
Options: 'FDR', 'p25', 'p50', 'p75' or any percentile in format p with number. To select for metric of comparison: -log10 FDR wilcoxon test corrected for multitesting with direction of ecdf shift, or difference from background for one of the percentiles from ecdf.

###### RegMode
Options: 'total', 'poly', 'translation'.  To select for mode of regulation, respectively: total mRNA log2FC, poly mRNA log2FC or translation log2 FC

###### pdfName
Optional, name of the output pdf file with plot

### 2. Features analysis:

In this part of analysis, ..... . In order to run all the following steps, anota2seq object and annotation file are required....

Run function retrieveFormatData to prepare annotation file. There are a number of options to choose from. It is possible to use existing pre-prepared files. At the moment only human and mouse are available (for other species use createFromFile option in function retrieveFormatData. To check existing versions run:

    checkAvailableVersions(species)


__Output:__ Prints out available versions for given species

###### species
Name of the species, for example human or mouse. 

It is possible to create new annotation files based on latest ncbi refseq release (option= 'create'). It is automatic process for human and mouse. For other species option = 'createFromFile' is required and files rna.fna.gz, rna.gbff.gz and genomic.gff.gz needs to be manually downloaded from ncbi refseq and provided. Once created or as an another option, it is possible simply load your own version in appropriate format, i.e. data.frame for transcripts IDs, gene IDs, and sequences for 5'UTR, CDS and 3'UTR with respective column names: id, geneID, UTR5_seq, CDS_seq, UTR3_seq.


    retrieveFormatData(source, species=NULL, version=NULL, customFile=NULL, rna_gbff_file=NULL, rna_gbff_file=NULL, rna_fa_file=NULL,genomic_gff_file=NULL)


__Output:__ Annotation file in appropriate format (loaded and written out if create or createFromFile), i.e. data.frame for transcripts IDs, gene IDs, and sequences for 5'UTR, CDS and 3'UTR with respective column names: id, geneID, UTR5_seq, CDS_seq, UTR3_seq

###### source
To choose from: load, custom, create or createFromFile

###### version
Loading available version of pre-existing annotation if source=='load'. Use function checkAvailableVersions(species) to check available releases

###### species
human and mouse only at the moment, it is required if source == 'create' or 'load'

###### customFile
path to the file, it is required if source == 'custom'. The file should be in format: data.frame for transcripts IDs, gene IDs, and sequences for 5'UTR, CDS and 3'UTR with respective column names: id, geneID, UTR5_seq, CDS_seq, UTR3_seq

###### rna_gbff_file
path to the file, it is required if source == 'createFromFile'. It can be downloaded from ncbi refseq.

###### rna_fa_file
path to the file, it is required if source == 'createFromFile'. It can be downloaded from ncbi refseq.

###### genomic_gff_file
path to the file, it is required if source == 'createFromFile'. It can be downloaded from ncbi refseq.

It is also possible to replace sequences if information from nanoCage or quantSeq exists. 

    adjustSeq(annot, region, adjObj, keepAll = FALSE, excl=FALSE)


__Output:__ Annotation file with adjusted sequenes in the same format as annotationfile (loaded and written out if create or createFromFile), i.e. data.frame for transcripts IDs, gene IDs, and sequences for 5'UTR, CDS and 3'UTR with respective column names: id, geneID, UTR5_seq, CDS_seq, UTR3_seq

###### annot
annotation file prepared by function: retrieveFormatData

###### region
UTR5 or UTR3

###### adjObj
data.frame object with refseq transcript ids and sequences (with column names: 'id', 'sequence')

###### keepAll
Default FALSE. For transcripts that have multiple isoforms, remove these for a given gene that do not have information about sequences in adjObj. The idea is to have only adjusted sequences for these genes when there is "better" information and not from database. Alternatively if TRUE it will also keep not adjusted isoforms for these genes alongside adjusted.

###### excl
Default FALSE, it keeps all genes/isoforms. If TRUE it will exclude genes/isoforms that are not in adjObj.

#### 2.1 Codon usage analysis

In this step, codon (or alternatively AA) usage is conducted. First function apply chi-square test and plots heatmap for all desired modes of regulation. Then, in addition to s additonal pairwise plots, for each codon odds ratio for each selected pair is calcuated and plotted against frequency of given codon (or AA). The codons (or AA) of interest are these with high odds ratio and high frequency. Applying second function (codonCalc) for these selected codons (or AA) it is possible to calculate overall counts/frequency for all genes that might be used in feature integration step. 

    codonUsage(ads, analysis, regulation, contrast, comparisons=NULL, type='sequence', annotType='ccds', source='load', species=NULL, annot=NULL, subregion=NULL,  subregionSel, selection='random', pAdj = 0.01, plotHeatmap=TRUE, pdfName=NULL)

__Output:__ Number of plots and codon/AA counts table for all the codons in all expressed genes 

###### ads
anota2seq object

###### analysis
'codon' or 'AA' to choose

###### regulation
desired names of regulation modes from anota2seq object to be compared. To choose from:  "translationUp", "translationDown", "bufferingmRNAUp", "bufferingmRNADown", "mRNAAbundanceUp", "mRNAAbundanceDown". It is possible to cross compare the same mode of regulation across different models within the same anota2seq object, by just providing twice the same name of regulation and stating correct contrasts. (for example c("translationUp","translationDown","translationUp")

###### contrast
respective contrasts in anota2seq object for each desired regulation mode (for example c(1,1,2) would indicate that first two modes of regulation provided above are from contrast 1 and third one from contrast 2)

###### comparisons
Indicate for which regulation pairs provide statistics by providing list of vectors (for example list(c(1,2),c(2,4))).

###### type
Default 'sequence'. There is an option to use ribosomal profiling by selection 'riboSeq'

###### annotType
Default 'ccds'. There is an option to use refSeq sequences. 

###### source
Default 'load'. it will load preprepared sequences for selected species. However, there is an option to automatically create new 'create' or load custom 'custom'

###### species
human and mouse only at the moment, it is required if source == 'create' or 'load'

###### annot
it is only required when annotType=='refSeq'. Annotation file prepared by function: retrieveFormatData

###### subregion
optional. number of nucleotides from start if positive or end if negative.

###### subregionSel
required if subregion is not null. 'select' or 'exclude'. whether subregion should be selected or excluded.

###### selection
Default: random, to choose from  'longest','shortest','random'. There are multiple isoforms for given gene, this decides whether the longest or the shortes for given region will be taken into consideration or simply random one.

###### pAdj
Default 0.01, threshold for chi-square test

###### plotHeatmap
Default TRUE, it plots heatmap of residuals

###### pdfName
Optional, name of the output pdf file with plot


    codonCalc(codonIn, analysis, featsel, unit='count', ads, regulation, contrast, comparisons=NULL, pdfName)

__Output:__ ecdf plots for comparisons and vector of counts/frequency for all expressed genes

###### codonIn
Output from codonUsage function

###### analysis
'codon' or 'AA' to choose

###### featsel
vector of desired codons/AA with high odds ratio and high frequency

###### unit
using 'count' or 'freq' for calculation

###### ads
anota2seq object

###### regulation
desired names of regulation modes from anota2seq object to be compared. To choose from:  "translationUp", "translationDown", "bufferingmRNAUp", "bufferingmRNADown", "mRNAAbundanceUp", "mRNAAbundanceDown". It is possible to cross compare the same mode of regulation across different models within the same anota2seq object, by just providing twice the same name of regulation and stating correct contrasts. (for example c("translationUp","translationDown","translationUp")

###### contrast
respective contrasts in anota2seq object for each desired regulation mode (for example c(1,1,2) would indicate that first two modes of regulation provided above are from contrast 1 and third one from contrast 2)

###### comparisons
Indicate for which regulation pairs provide statistics by providing list of vectors (for example list(c(1,2),c(2,4))).

###### pdfName
Optional, name of the output pdf file with plot

#### 2.2 Length analysis

In this step, lengths of desired region of mRNA is compared.

    lengthAnalysis(ads, annot, regulation, contrast, comparisons=NULL, region, selection='random', plotType='boxplot', pdfName=NULL)

__Output:__ Desired plot and vector of lengths for all expressed genes

###### ads
anota2seq object

###### annot
annotation file prepared by function: retrieveFormatData

###### regulation
desired names of regulation modes from anota2seq object to be compared. To choose from:  "translationUp", "translationDown", "bufferingmRNAUp", "bufferingmRNADown", "mRNAAbundanceUp", "mRNAAbundanceDown". It is possible to cross compare the same mode of regulation across different models within the same anota2seq object, by just providing twice the same name of regulation and stating correct contrasts. (for example c("translationUp","translationDown","translationUp")

###### contrast
respective contrasts in anota2seq object for each desired regulation mode (for example c(1,1,2) would indicate that first two modes of regulation provided above are from contrast 1 and third one from contrast 2)

###### comparisons
Optional. Indicate for which regulation pairs provide statistics by providing list of vectors (for example list(c(0,1),c(0,2), c(1,2)). 0 is background. It will compare regulation 1 to background and regulation 2 to background and regulation 1 and regulation 2 between each other).

###### region
To choose from 3'UTR, CDS, 5'UTR

###### selection
Default: random, to choose from  'longest','shortest','random'. There are multiple isoforms for given gene, this decides whether the longest or the shortes for given region will be taken into consideration or simply random one.

###### plotType
Default: boxplot, to choose from 'boxplot', 'violin', 'ecdf'

###### pdfName
Optional, name of the output pdf file with plot

#### 2.3 Content analysis

In this step, nucleotide content of choice in desired region of mRNA is compared.

    contentAnalysis(ads, annot, content, regulation, contrast, comparisons=NULL, region, subregion=NULL, subregionSel, selection='random', plotType='boxplot', pdfName=NULL)

__Output:__ Desired plot and vector of %content for all expressed genes

###### ads
anota2seq object

###### annot
annotation file prepared by function: retrieveFormatData

###### content
to select from  A,T,G,C in any combination

###### regulation
desired names of regulation modes from anota2seq object to be compared. To choose from:  "translationUp", "translationDown", "bufferingmRNAUp", "bufferingmRNADown", "mRNAAbundanceUp", "mRNAAbundanceDown". It is possible to cross compare the same mode of regulation across different models within the same anota2seq object, by just providing twice the same name of regulation and stating correct contrasts. (for example c("translationUp","translationDown","translationUp")

###### contrast
respective contrasts in anota2seq object for each desired regulation mode (for example c(1,1,2) would indicate that first two modes of regulation provided above are from contrast 1 and third one from contrast 2)

###### comparisons
optional. Indicate for which regulation pairs provide statistics by providing list of vectors (for example list(c(0,1),c(0,2), c(1,2)). 0 is background. It will compare regulation 1 to background and regulation 2 to background and regulation 1 and regulation 2 between each other).

###### region
to choose from 3'UTR, CDS, 5'UTR

###### subregion
optional. number of nucleotides from start if positive or end if negative.

###### subregionSel
required if subregion is not null. 'select' or 'exclude'. whether subregion should be selected or excluded.

###### selection
Default: random, to choose from  'longest','shortest','random'. There are multiple isoforms for given gene, this decides whether the longest or the shortes for given region will be taken into consideration or simply random one.

###### plotType
Default: boxplot, to choose from 'boxplot', 'violin', 'ecdf'

###### pdfName
Optional, name of the output pdf file with plot

#### 2.4 Discovery of motifs

In this step, motif discovery using streme and matching to desired known motifs in database of choice.

    motifAnalysis(ads, annot, regulation, contrast, geneList=NULL, region, subregion=NULL, subregionSel, selection='random', stremeThreshold = 0.05, stremeName=NULL, tomtomName=NULL, tomtom_database)

__Output:__ Dataframe with streme and tomotom results. Fasta files with sequences for control and regulated genes is also written out. In order to plot desired motif: motifOut[select here row numbers like normal indexing ,] %>% to_list() %>% view_motifs()

###### ads
anota2seq object

###### annot
annotation file prepared by function: retrieveFormatData

###### regulation
desired name of regulation mode from anota2seq object to be compared (select only one mode). To choose from:  "translationUp", "translationDown", "bufferingmRNAUp", "bufferingmRNADown", "mRNAAbundanceUp", "mRNAAbundanceDown".

###### contrast
select contrast in anota2seq object for desired regulation mode (select only one contrast)

###### region
to choose from 3'UTR, CDS, 5'UTR

###### subregion
optional. number of nucleotides from start if positive or end if negative.

###### subregionSel
required if subregion is not null. 'select' or 'exclude'. whether subregion should be selected or excluded.

###### selection
Default: random, to choose from  'longest','shortest','random'. There are multiple isoforms for given gene, this decides whether the longest or the shortes for given region will be taken into consideration or simply random one.

###### stremeThreshold
Default 0.05. I think it is max possible.

###### stremeName
optional. name for output folder for streme analysis

###### tomtomName
optional. name for output folder for tomtom analysis

###### tomtom_database
path to .meme file containing database (including name of the file). All databases can be downloaded from here: https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.21.tgz

#### 2.5 Presence of motifs

In this step, quantitative presence of motifs in the sequences is calculated. There is also an option to predict presence of g-quadruplexes

    contentMotifs(ads, annot, regulation, contrast, comparisons=NULL, motif, len=1, min_score=47, region, subregion=NULL, subregionSel, selection='random',pdfName=NULL)

__Output:__  Plot and vector of %content for all expressed genes

###### ads
anota2seq object

###### annot
annotation file prepared by function: retrieveFormatData

###### regulation
desired name of regulation mode from anota2seq object to be compared. To choose from:  "translationUp", "translationDown", "bufferingmRNAUp", "bufferingmRNADown", "mRNAAbundanceUp", "mRNAAbundanceDown". It is possible to cross compare the same mode of regulation across different models within the same anota2seq object, by just providing twice the same name of regulation and stating correct contrasts. (for example c("translationUp","translationDown","translationUp")

###### contrast
respective contrasts in anota2seq object for each desired regulation mode (for example c(1,1,2) would indicate that first two modes of regulation provided above are from contrast 1 and third one from contrast 2)

###### comparisons
Indicate for which regulation pairs provide statistics by providing list of vectors (for example list(c(1,2),c(2,4))).

###### motif
motif sequence. it can be IUPAC code or by using [ ] for ambiguities. Use 'G4' for g-quadruplexes. 

###### len
min distance between motifs. not required if 'G4'

###### min_score
default 47. quality of the prediction, required for 'G4'. Read more about: https://bioconductor.org/packages/release/bioc/html/pqsfinder.html

###### region
to choose from 3'UTR, CDS, 5'UTR

###### subregion
optional. number of nucleotides from start if positive or end if negative.

###### subregionSel
required if subregion is not null. 'select' or 'exclude'. whether subregion should be selected or excluded.

###### selection
Default: random, to choose from  'longest','shortest','random'. There are multiple isoforms for given gene, this decides whether the longest or the shortes for given region will be taken into consideration or simply random one.

###### pdfName
Optional, name of the output pdf file with plot

#### 2.6 Signature

In this step, signature as input for feature integration is created (vector of  absence/presence in signature of genes)

    signCalc(ads, signature)

__Output:__  Vector for all expressed genes (0 or 1) absent/present in signature

###### ads
anota2seq object

###### signature
vector of gene names in given signature

#### 2.7 Folding energy

Here, differences in folding energy between modes of regulation can be assessed. As the process of calculating folding energy is long, it is possible to use the function only to preprepare the input file from fasta or based on annot (it requires mfold). It is also possible just to use already prepared files (for mouse and human).


    foldEnergyAnalysis(source='load', version, species, region=NULL,fromFasta=FALSE,customFile, onlyRun=FALSE, ads, regulation, contrast, comparisons, annot, selection, plotType='boxplot', resid=TRUE, pdfName=NULL)

__Output:__  Plots and vector of folding energies

###### source
Default 'load'. it will load preprepared fold energies for selected species. However, there is an option to automatically create new 'create' based on selected region in annot or fromFasta ( fasta files with sequences and transcript ids as names) or load custom 'custom' file with preprepared fold energies calculations in format: transcirpt ID, folding energy, length of sequence (with column names: 'id','fold_energy','length').

###### version
Loading available version of pre-existing annotation if source=='load'. Use function checkAvailableVersions(species) to check available releases

###### species
human and mouse only at the moment, it is required for source== 'load'

###### region
to choose from 3'UTR, CDS, 5'UTR. It is not required for custom options

###### fromFasta
Default is FALSE, it is an option to run mfold on custom fasta file.

###### customFile
path and file. It is required for source=='custom' or if source=='create' and fromFasta==TRUE. In former case it must be the file with calculated folding energy in format: transcirpt ID, folding energy, length of sequence (with column names: 'id','fold_energy','length') or in latter fasta files with sequences and transcript ids as names

###### onlyRun
Default FALSE, it will load or create fold energy input file using mfold and run analysis comparing differnt modes of regulation. There is an option to just run mfold to prepare fold energy file without analysis.

###### ads
anota2seq object

###### regulation
desired name of regulation mode from anota2seq object to be compared (select only one mode). To choose from:  "translationUp", "translationDown", "bufferingmRNAUp", "bufferingmRNADown", "mRNAAbundanceUp", "mRNAAbundanceDown".

###### contrast
respective contrasts in anota2seq object for each desired regulation mode (for example c(1,1,2) would indicate that first two modes of regulation provided above are from contrast 1 and third one from contrast 2)

###### comparisons
Indicate for which regulation pairs provide statistics by providing list of vectors (for example list(c(1,2),c(2,4))).

###### annot
annotation file prepared by function: retrieveFormatData

###### selection
Default: random, to choose from  'longest','shortest','random'. There are multiple isoforms for given gene, this decides whether the longest or the shortes for given region will be taken into consideration or simply random one.

###### plotType
Default: boxplot, to choose from 'boxplot', 'ecdf'

###### resid
Default TRUE, residuals after correction for length or option to choose raw fold energy values

###### pdfName
Optional, name of the output pdf file with plot

#### 2.8 miRNA enrichment

#### 2.9 uORFs presence

In this step differences in presence of uORFs are calculated. There is an option to look for uORFs only in 5'UTR (completely contained in) or starting in 5'UTR but with possiblity of having stop codon in CDS or 3'UTR. There is an option to use different start codons in differnt Kozak contexts.

    uorf_analysis(ads, regulation, contrast, comparisons=NULL, annot, selection, onlyUTR5=FALSE, startCodon='ATG', KozakContext='strong', pdfName=NULL)

__Output:__  Plot and vector with number of uORFs per gene
 
###### ads
anota2seq object

###### regulation
desired name of regulation mode from anota2seq object to be compared (select only one mode). To choose from:  "translationUp", "translationDown", "bufferingmRNAUp", "bufferingmRNADown", "mRNAAbundanceUp", "mRNAAbundanceDown".

###### contrast
respective contrasts in anota2seq object for each desired regulation mode (for example c(1,1,2) would indicate that first two modes of regulation provided above are from contrast 1 and third one from contrast 2)

###### comparisons
Indicate for which regulation pairs provide statistics by providing list of vectors (for example list(c(1,2),c(2,4))).

###### annot
annotation file prepared by function: retrieveFormatData

###### selection
Default: random, to choose from  'longest','shortest','random'. There are multiple isoforms for given gene, this decides whether the longest or the shortes for given region will be taken into consideration or simply random one.

###### onlyUTR5
Default FALSE. Looking for uORFs starting in 5'UTR but with possiblity of having stop codon in CDS or 3'UTR. If TRUE , only in 5'UTR (completely contained in).

###### startCodon
Default 'ATG, but it is possible to provide any start codon

###### KozakContext
Default is: strong (i.e [AG]..startCodonG). But there are options to choose 'adequate1' ([AG]..startCodon.) or 'adequate2' (...startCodonG) or 'weak' (...startCodon.) and finally 'any' uORFs in any kozak context

###### pdfName
Optional, name of the output pdf file with plot

#### 2.10 Feature integration 

All steps in Feature Analysis provided in output calculations for all genes that can be used as input for feature integration that will calculate importance of each.


    featureIntegration(ads, contrast, RegMode='translation' , pdfName=NULL)

__Output:__  plot table with all the calculations

###### ads
anota2seq object

###### contrast
select contrast in anota2seq object (select only one contrast)

###### RegMode
Default='translation'. select mode of regulation, to choose from: translation, buffering, translatedmRNA, totalmRNA

###### features
list of features for analysis given in format of named list

###### pdfName
Optional, name of the output pdf file with plot



------------------------------------------------------------------------

## Contacts

krzysztof.szkop@gmail.com

