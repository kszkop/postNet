\name{signCalc}
\alias{signCalc}
\title{Create gene signatures for use in feature integration analysis
}

\description{The \code{signCalc} function converts lists of gene IDs into signatures that are compatible as inputs into downstream
\link[postNet]{featureIntegration} models.
}

\usage{
signCalc(ptn, 
        signatures)
}

\arguments{
  \item{ptn}{A \code{postNetData} object.
  }
  \item{signatures}{A named list of character vectors containing gene IDs to be converted into signatures
  that can be used as input with the \link[postNet]{featureIntegration} function.
  }
}

\value{Named list of vectors indicating inclusion in different signatures for all genes in the \code{postNetData} object. 
Genes in the signature that are included in the \code{postNetData} object are encoded as 1, and those that are absent are 
encoded as a 0. 
}

\seealso{
\link[postNet]{featureIntegration}
}

\examples{
# load example data:
 data(postNetExample)
    
# Initialize the postNetData object with custom gene lists:
 translationUp <- postNetExample$geneList[[1]]
 translationDown <- postNetExample$geneList[[2]]
 myGenes <- list('translationUp' = translationUp, 'translationDown' = translationDown)
 myBg <- postNetExample$background
 myEffect <- postNetExample$effect

 ptn <- postNetStart(
   geneList = myGenes,
   geneListcolours = c("#7FB7BE","#DB7F67"),
   customBg = myBg,
   effectMeasure = myEffect,
   selection = "random",
   setSeed = 123,
   source = "load", 
   species = "human"
 )
 
# Create signatures of interest to be used in featureIntegration models
 addSign < -list()
 addSign[["TOP_mRNAs"]] <- c("RPL4", "RPL5","EEF1A1","EIF3D")
 addSign[["Hypoxia_Induced"]] <- c("VEGFA","LOX","ENO1","PDK1","PGK1","HIF1")

 mySign <- signCalc(ptn = ptn, signatures = addSign)
 str(mySign)
}
