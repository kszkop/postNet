\name{mouseSignatures}
\alias{mouseSignatures}
\docType{data}
\title{Mouse gene signatures of translational control
}
\description{The \code{mouseSignatures} dataset is a list of 13 gene signatures derived from published studies profiling mRNA
translation downstream of various pathways and factors, or mRNAs containing features such as 5'TOP motifs. The list is not intended to be exhaustive,
but provides examples of how gene signatures can be used in a \pkg{postNet} analysis to gain insights into post-transcriptional regulation
of gene expression.}

\usage{data("mouseSignatures")}
\format{
  The format is:
List of 13
 $ Lee_etal_2014_senescence_Up        : chr [1:598] "Zfy1" "Zfy2" "Asxl3" "Zfp397" ...
 $ Lee_etal_2014_senescence_Down      : chr [1:714] "Cdh2" "Zfp185" "Cenpm" "Has2" ...
 $ Larsson_etal_2007_high4E           : chr [1:1251] "Asl" "Emd" "Fah" "Galt" ...
 $ Karampelias_etal_2022_p4E_transUp  : chr [1:252] "Prkaa1" "Uprt" "Ttc33" "Lratd2" ...
 $ Karampelias_etal_2022_p4E_transDown: chr [1:57] "Urb1" "Odc1" "Ccdc85b" "Sbf1" ...
 $ Gandin_etal_2016_mTOR_transUp      : chr [1:1487] "Zfy1" "Zfy2" "Isca2" "Rps14" ...
 $ Gandin_etal_2016_mTOR_transDown    : chr [1:1765] "Lpin1" "Galnt1" "Atrx" "Tbc1d31" ...
 $ Cockman_etal_2020_classicTOP       : chr [1:70] "Rpl37" "Rps14" "Hnrnpa1" "Rpl39" ...
 $ Han_etal_2013_ATF4                 : chr [1:254] "Akap2" "Apbb2" "Apobec1" "Slc7a1" ...
 $ Guan_etal_2017_Tg1_transUp         : chr [1:1336] "Mcm7" "Bcl2" "Usp33" "Gm13212" ...
 $ Guan_etal_2017_Tg1_transDown       : chr [1:774] "Mapre3" "Mettl24" "Isg15" "Idh3g" ...
 $ Baird_etal_2014_Tg_transUp         : chr [1:205] "Atf5" "Ibtk" "Sf3b1" "Cdk5rap2" ...
 $ Baird_etal_2014_Tg_transDown       : chr [1:230] "Gas1" "Angptl2" "Ppp3cb" "Scara3" ...
}
\details{The \code{mouseSignatures} dataset contains a curated list of gene signatures relevent to understanding possible pathways and
factors involved in regulation of mRNA translation, with a specific focus on mTOR and integrated stress response (ISR)-sensitive translation. 
These signatures have been compiled either directly from published studies, or from analyses of published datasets (sources and references 
described in detail below). Note that if the orignal source of the gene signature was from human, these gene IDs have been converted to mouse. 

These gene signatures can be used with the \link[postNet]{plotSignatures}, \link[postNet]{plotSignatures_ads}, 
and \link[postNet]{signaturesHeatmap} functions to perform statistical analyses and visualizations of how each
gene signature is regulated in a given dataset. Furthermore, gene signatures can also be provided as features 
in the \link[postNet]{featureIntegration} function which performs step-wise regression modelling or random forest
classification to identify features that are most important in explaining post-transcriptional regulation in a given
dataset. 
}

\source{
\strong{High eIF4E-dependent translation:} \cr
 \itemize{
 \item \strong{Signature Information}:{
    \itemize{
        \item Signature Name: Larsson_etal_2007_high4E
        \item Source: HMEC cells immortalized with hTERT
        \item Species: Human
        \item Conditions: HMEC/hTERT + eIF4E over-expression (OE), HMEC/hTERT
        \item Replicates: 2 Reps/Condition (HMEC/hTERT + eIF4E OE polyome-associated mRNA had 3 samples)
        \item Comparison: HMEC/hTERT + eIF4E OE vs. HMEC/hTERT 
  }
  }
  \item \strong{Publication Information}:{
    \itemize{
        \item Title: Eukaryotic translation initiation factor 4E induced progression of primary human mammary epithelial cells along the cancer pathway is associated with targeted translational deregulation of oncogenic drivers and inhibitors
        \item Method: Polysome fractionation followed by microarray 
        \item Data Source: Data obtained from supplementary materials (see: Table S2 of publication)
        \item Analysis Method: Significance analysis of microarrays algorithm, genes with FDR < 0.10.
        \item First Author,Last Author, Year:  Larsson, Polunovsky, 2007
        \item PMID: 17638893
        \item GEO Accession: GSE6043
        \item DOI: 10.1158/0008-5472.CAN-07-0752
    }
    }
 \item \strong{Dataset Summary}:{
 Genes classified as translationally activated with high eIF4E expression in HMEC/hTERT + eIF4E OE compared to HMEC/hTERT.  
  }
}
  
\strong{Genes translationally regulated downstream of high phospho-eIF4E:} \cr
 \itemize{
 \item \strong{Signature Information}:{
    \itemize{
        \item Signature Name(s): Karampelias_etal_2022_p4E_transUp, Karampelias_etal_2022_p4E_transDown
        \item Source: PANC1 cells
        \item Species: Human
        \item Conditions: Cercosporamide (5 μM), CID661578 (40 μM), DMSO: 6 hrs
        \item Replicates: 4 Reps/Condition
        \item Comparison: Cercosporamide vs. DMSO, CID661578 vs. DMSO
  }
  }
  \item \strong{Publication Information}:{
    \itemize{
        \item Title: MNK2 deficiency potentiates β-cell regeneration via translational regulation
        \item Method: Polysome fractionation followed by RNA sequencing 
        \item Data Source: Data obtained from supplementary materials (see: Supplementary Data 3 of publication)
        \item Analysis Method: Anota2seq algorithm (version 1.14.0). The following thresholds were applied within the anota2seqRun function: 
        maxPAdj = 0.15; deltaP = log2(1); deltaT = log2(1); deltaPT = log2(1.2); deltaTP = log2(1.2); maxSlopeTranslation = 2; 
        minSlopeTranslation = −1; minSlopeBuffering = −2; maxSlopeBuffering = 1. Replicate was included in the anota2seq model using the ‘batchVec’ 
        parameter to account for batch effects.
        \item First Author, Last Author, Year:  Karampelias, Andersson, 2022
        \item PMID: 35697798
        \item GEO Accession: GSE200477
        \item DOI: 10.1038/s41589-022-01047-x
    }
    }
 \item \strong{Dataset Summary}:{
 The intersect of genes classified as translationally activated (transUp) or suppressed (transDown) following cercosporamide and CID661578 treatment
 compared to DMSO control in PANC1 cells. Both cercosporamide and CID661578 induce high phosphorylation of eIF4E in PANC1 cells.
  }
}

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

Larsson O, Li S, Issaenko OA, Avdulov S, Peterson M, Smith K, Bitterman PB, Polunovsky VA. Eukaryotic translation initiation factor 4E induced 
progression of primary human mammary epithelial cells along the cancer pathway is associated with targeted translational deregulation of oncogenic 
drivers and inhibitors. Cancer Res. 2007 Jul 15;67(14):6814-24. doi: 10.1158/0008-5472.CAN-07-0752. PMID: 17638893.

Karampelias C, Watt K, Mattsson CL, Ruiz ÁF, Rezanejad H, Mi J, Liu X, Chu L, Locasale JW, Korbutt GS, Rovira M, Larsson O, Andersson O. 
MNK2 deficiency potentiates β-cell regeneration via translational regulation. Nat Chem Biol. 2022 Sep;18(9):942-953. doi: 10.1038/s41589-022-01047-x. 
Epub 2022 Jun 13. PMID: 35697798; PMCID: PMC7613404.

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
 data(mouseSignatures)
 str(mouseSignatures)

}
\keyword{datasets}
