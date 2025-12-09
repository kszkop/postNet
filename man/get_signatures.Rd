\name{get_signatures}
\alias{get_signatures}
\title{Retrieve gene signatures provided in the postNetParcel}

\description{The \code{get_signatures} function retrieves gene signatures of interest provided with the package for the indicated species. 
}

\usage{
get_signatures(species)
}

\arguments{
  \item{species}{A string to select the species of interest. Currently, \code{"human"} and \code{"mouse"} are supported.
}
}

\details{The data retrieved contains a list of gene signatures relevant to understanding possible pathways and mRNA features involved in regulation of mRNA translation, with a specific focus on mTOR and integrated stress response (ISR)-sensitive translation. These signatures have been compiled either directly from published studies, or from analyses of published datasets (sources and references described in detail below). Note that if the original source of the gene signature was from human, these gene IDs have been converted to mouse, and vice versa. 

These gene signatures can be used with the \link[postNet]{plotSignatures}, \link[postNet]{plotSignatures_ads}, and \link[postNet]{signaturesHeatmap} functions to perform statistical analyses and visualizations of how each gene signature is regulated in a given dataset. Furthermore, gene signatures can also be provided as features in the \link[postNet]{featureIntegration} function which performs forward stepwise regression modelling or Random Forest classification to identify features that are most important in explaining post-transcriptional regulation in a given dataset. The \link[postNet]{signCalc} function can be used to convert gene signatures into a format compatible with \link[postNet]{featureIntegration}.   
}

\value{
The format is: \cr
List of 5 \cr
 $ Gandin_etal_2016_mTOR_transUp      : chr [1:1487] "Zfy1" "Zfy2" "Isca2" "Rps14" ...
 $ Gandin_etal_2016_mTOR_transDown    : chr [1:1765] "Lpin1" "Galnt1" "Atrx" "Tbc1d31" ...
 $ Cockman_etal_2020_classicTOP       : chr [1:70] "Rpl37" "Rps14" "Hnrnpa1" "Rpl39" ...
 $ Guan_etal_2017_Tg1_transUp         : chr [1:1336] "Mcm7" "Bcl2" "Usp33" "Gm13212" ...
 $ Guan_etal_2017_Tg1_transDown       : chr [1:774] "Mapre3" "Mettl24" "Isg15" "Idh3g" ...
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
        \item Method: Polysome fractionation followed by RNA sequencing 
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
If you use these gene signatures in your analysis, please cite: \cr

Gandin V, Masvidal L, Cargnello M, Gyenis L, McLaughlan S, Cai Y, Tenkerian C, Morita M, Balanathan P, Jean-Jean O, Stambolic V, Trost M, Furic L, 
Larose L, Koromilas AE, Asano K, Litchfield D, Larsson O, Topisirovic I. mTORC1 and CK2 coordinate ternary and eIF4F complex assembly. Nat Commun. 
2016 Apr 4;7:11127. doi: 10.1038/ncomms11127. PMID: 27040916; PMCID: PMC4822005.

Cockman E, Anderson P, Ivanov P. TOP mRNPs: Molecular Mechanisms and Principles of Regulation. Biomolecules. 2020 Jun 27;10(7):969. 
doi: 10.3390/biom10070969. PMID: 32605040; PMCID: PMC7407576.

Guan BJ, van Hoef V, Jobava R, Elroy-Stein O, Valasek LS, Cargnello M, Gao XH, Krokowski D, Merrick WC, Kimball SR, Komar AA, Koromilas AE, 
Wynshaw-Boris A, Topisirovic I, Larsson O, Hatzoglou M. A Unique ISR Program Determines Cellular Responses to Chronic Stress. Mol Cell. 
2017 Dec 7;68(5):885-900.e6. doi: 10.1016/j.molcel.2017.11.007. PMID: 29220654; PMCID: PMC5730339.
}

\seealso{
\link[postNet]{plotSignatures} \cr
\link[postNet]{plotSignatures_ads} \cr
\link[postNet]{plotSignaturesHeatmap} \cr
\link[postNet]{signCalc}
}

\examples{
# Load the human gene signatures: 
 humanSignatures <- get_signatures("human")
 str(humanSignatures)

}
\keyword{datasets}

