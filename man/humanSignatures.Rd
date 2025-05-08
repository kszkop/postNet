\name{humanSignatures}
\alias{humanSignatures}
\docType{data}
\title{Human gene signatures of translational control
}

\description{The \code{humanSignatures} dataset is a list of example gene signatures derived from published studies profiling mRNA
translation downstream of mTOR signalling and activation of the integrated stress response (ISR). In addition, an example of a signature comprised of mRNAs containing 
5'TOP motifs is also included. The list is not intended to be exhaustive, but provides examples of how gene signatures can be used in a \pkg{postNet} analysis to gain 
insights into post-transcriptional regulation of gene expression.} 

\usage{data("humanSignatures")}
\format{
  The format is:
List of 13
 $ Lee_etal_2014_senescence_Up        : chr [1:632] "TXNIP" "OLR1" "PAPPA" "KIAA1199" ...
 $ Lee_etal_2014_senescence_Down      : chr [1:723] "SFPQ" "MRE11A" "C5" "PDXP" ...
 $ Larsson_etal_2007_high4E           : chr [1:2146] "EPHB2" "A1BG" "AAAS" "AAMP" ...
 $ Karampelias_etal_2022_p4E_transUp  : chr [1:270] "HK2" "MAP3K20-AS1" "CXCL8" "ENPP4" ...
 $ Karampelias_etal_2022_p4E_transDown: chr [1:99] "COX3" "CYTB" "ATP6" "COX1" ...
 $ Gandin_etal_2016_mTOR_transUp      : chr [1:1569] "ZNF467" "PHLDA2" "PYCARD" "FHOD1" ...
 $ Gandin_etal_2016_mTOR_transDown    : chr [1:1721] "DST" "OTULINL" "MME" "BRD8" ...
 $ Cockman_etal_2020_classicTOP       : chr [1:92] "RPSA" "RPS2" "RPS3" "RPS3A" ...
 $ Han_etal_2013_ATF4                 : chr [1:225] "USP11" "CRLS1" "APBB2" "COL18A1" ...
 $ Guan_etal_2017_Tg1_transUp         : chr [1:1192] "AR" "ATP7A" "BRCA2" "LYST" ...
 $ Guan_etal_2017_Tg1_transDown       : chr [1:731] "AGA" "BCKDHB" "CST3" "CYBA" ...
 $ Baird_etal_2014_Tg_transUp         : chr [1:179] "HADHA" "NPC1" "GPD2" "ATRX" ...
 $ Baird_etal_2014_Tg_transDown       : chr [1:194] "HMBS" "SGCB" "SMN1" "SMN2" ...
}

\details{The \code{humanSignatures} dataset contains a list of gene signatures relevent to understanding possible pathways and
mRNA features involved in regulation of mRNA translation, with a specific focus on mTOR and integrated stress response (ISR)-sensitive translation. 
These signatures have been compiled either directly from published studies, or from analyses of published datasets (sources and references 
described in detail below). Note that if the orignal source of the gene signature was from mouse, these gene IDs have been converted to human. 

These gene signatures can be used with the \link[postNet]{plotSignatures}, \link[postNet]{plotSignatures_ads}, 
and \link[postNet]{signaturesHeatmap} functions to perform statistical analyses and visualizations of how each
gene signature is regulated in a given dataset. Furthermore, gene signatures can also be provided as features 
in the \link[postNet]{featureIntegration} function which performs step-wise regression modelling or random forest
classification to identify features that are most important in explaining post-transcriptional regulation in a given
dataset. 
}

\source{
\strong{Genes translationally regulated downstream of mTOR:} \cr
 \itemize{
 \item \strong{Signature Information}:{
    \itemize{
        \item Signature Name(s): Gandin_etal_2016_mTOR_transUp, Gandin_etal_2016_mTOR_transDown
        \item Source: MCF7 cells
        \item Species: Human
        \item Conditions: Insulin (4.2 nM) + torin1 (250 nM), Insulin (4.2 nM) + DMSO; 4 hrs
        \item Replicates: 4 Reps/Condition 
        \item Comparison: Insulin + torin1 vs. Insulin
  }
  }
  \item \strong{Publication Information}:{
    \itemize{
        \item Title: mTORC1 and CK2 coordinate ternary and eIF4F complex assembly  
        \item Method: Polysome fractionation followed by microarray 
        \item Data Source: Anota2seq analysis of microarray data deposited in GSE76766.   
        \item Analysis Method: Anota2seq algorithm (version 1.14.0). The following thresholds were applied within the anota2seqRun function: 
        maxPAdj = 0.15; deltaP = log2(1.2); deltaT = log2(1.2); deltaPT = log2(1.2); deltaTP = log2(1.2); maxSlopeTranslation = 2; 
        minSlopeTranslation = −1; minSlopeBuffering = −2; maxSlopeBuffering = 1. Replicate was included in the anota2seq model using the ‘batchVec’ 
        parameter to account for batch effects. 
        \item First Author, Last Author, Year:  Gandin, Topisirovic, 2016
        \item PMID: 27040916
        \item GEO Accession: GSE76766
        \item DOI: 10.1038/ncomms11127
    }
    }
 \item \strong{Dataset Summary}:{
 Genes classified as translationally activated (transUp) or suppressed (transDown) upon inhibition of mTOR via torin1 treatment in insulin-stimulated MCF7 cells.
  }
}

\strong{mRNAs containing 5'terminal oligopyrimidine (5'TOP) motifs:} \cr
 \itemize{
 \item \strong{Signature Information}:{
    \itemize{
        \item Signature Name: Cockman_etal_2020_classicTOP
        \item Source: N/A
        \item Species: Human
        \item Conditions: N/A
        \item Replicates: N/A
        \item Comparison: N/A
  }
  }
  \item \strong{Publication Information}:{
    \itemize{
        \item Title: TOP mRNPs: Molecular Mechanisms and Principles of Regulation
        \item Method: N/A 
        \item Data Source: Data obtained from supplementary materials (see: Appendix A Table 1A of publication)
        \item Analysis Method:
        \item First Author, Last Author, Year:  Cockman, Ivanov, 2020
        \item PMID: 32605040
        \item GEO Accession: N/A
        \item DOI: 10.3390/biom10070969
    }
    }
 \item \strong{Dataset Summary}:{
 A curated list of well-validated 5'TOP motif-containing mRNAs.
  }
}

\strong{Genes translationally regulated downstream of ISR activation:} \cr
 \itemize{
 \item \strong{Signature Information}:{
    \itemize{
        \item Signature Name(s): Guan_etal_2017_Tg1_transUp, Guan_etal_2017_Tg1_transDown
        \item Source: Mouse embryonic fibroblasts (MEFs)
        \item Species: Mouse
        \item Conditions: Thapsigargin (400 nM), Vehicle control (DMSO); 1 hr
        \item Replicates: 4 Reps/Condition 
        \item Comparison: Thapsigargin vs. DMSO
  }
  }
  \item \strong{Publication Information}:{
    \itemize{
        \item Title: A Unique ISR Program Determines Cellular Responses to Chronic Stress
        \item Method: Polysome fractionation followed by RNA seqeuncing 
        \item Data Source: 
        \item Analysis Method:
        \item First Author, Last Author, Year:  Guan, Hatzoglou, 2017
        \item PMID: 29220654
        \item GEO Accession: GSE90070
        \item DOI: 10.1016/j.molcel.2017.11.007
    }
    }
 \item \strong{Dataset Summary}:{
 Genes translationally activated (transUp) or suppressed (transDown) following activation of the ISR via thapsigargin treatment in MEFs. 
  }
}
}

\references{
If you use these gene signatures in you analysis, please cite: \cr

Gandin V, Masvidal L, Cargnello M, Gyenis L, McLaughlan S, Cai Y, Tenkerian C, Morita M, Balanathan P, Jean-Jean O, Stambolic V, Trost M, Furic L, 
Larose L, Koromilas AE, Asano K, Litchfield D, Larsson O, Topisirovic I. mTORC1 and CK2 coordinate ternary and eIF4F complex assembly. Nat Commun. 
2016 Apr 4;7:11127. doi: 10.1038/ncomms11127. PMID: 27040916; PMCID: PMC4822005.

Cockman E, Anderson P, Ivanov P. TOP mRNPs: Molecular Mechanisms and Principles of Regulation. Biomolecules. 2020 Jun 27;10(7):969. 
doi: 10.3390/biom10070969. PMID: 32605040; PMCID: PMC7407576.

Guan BJ, van Hoef V, Jobava R, Elroy-Stein O, Valasek LS, Cargnello M, Gao XH, Krokowski D, Merrick WC, Kimball SR, Komar AA, Koromilas AE, 
Wynshaw-Boris A, Topisirovic I, Larsson O, Hatzoglou M. A Unique ISR Program Determines Cellular Responses to Chronic Stress. Mol Cell. 
2017 Dec 7;68(5):885-900.e6. doi: 10.1016/j.molcel.2017.11.007. PMID: 29220654; PMCID: PMC5730339.
}

\examples{
# Load the mouse signatures 
 data(humanSignatures)
 str(humanSignatures)}

\keyword{datasets}
