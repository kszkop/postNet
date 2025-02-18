\docType{methods}
\name{ptn_selection}
\alias{ptn_selection,postNetData-method}
\title{Accessor for the isoform selection parameter from the 'selection' slot of a \code{postNetData} object
}

\description{The 'selection' slot holds the method chosen to select which mRNA isoforms are used in analyses. These can be one
  of: \code{'longest'}, \code{'shortest'}, or \code{'random'}.
}

\usage{
\S4method{ptn_selection}{postNetData}(x)
}

\arguments{
  \item{x}{a \code{postNetData} object.
  }
}

\value{A character specifying the method of selecting mRNA isoforms. Can be either \code{'longest'}, \code{'shortest'}, 
or \code{'random'})
}

\examples{

# load and create example data:
 data("postNetData")

 myGenes <- ptn_getExampleGeneList(ads)
 myBg <- ptn_getExampleBg(ads)
 myEffect <- ptn_getExampleEffect(ads)

# Initialize the postNetData object with custom gene lists:
 ptn <- postNetStart(
   geneList = myGenes,
   geneListcolours = c("#7FB7BE","#DB7F67"),
   customBg = myBg,
   effectMeasure = myEffect,
   selection = "random",
   #setSeed = 123
   source = "load", 
   species = "human"
 )

# Access mRNA isoform selection method
  ptn_selection(ptn)

}

\seealso{
\code{\link[postNet]{postNetStart}} for details on the different mRNA isoform selection options.
}

