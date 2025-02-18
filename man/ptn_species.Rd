\docType{methods}
\name{ptn_species}
\alias{ptn_species,postNetData-method}
\title{Accessor for the "species" slot of a \code{postNetData} object
}

\description{The "species" slot holds the source species of the RefSeq sequence annotations. Currently \code{"human"} and 
\code{"mouse"} are available if \code{source} is \code{"load"}, \code{"create"}, \code{"createFromSourceFiles"}, or 
\code{"createFromFasta"}. 
}

\usage{
\S4method{ptn_species}{postNetData}(x)
}

\arguments{
  \item{x}{a \code{postNetData} object.
  }
}

\value{A character specifying the species of origin for the RefSeq annotations used to generate the reference sequences.
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

# Access species of origin for the RefSeq annotation used to create the reference sequences
  ptn_species(ptn)

}

\seealso{
\code{\link[postNet]{postNetStart}} for details on available species, and providing custom annotations.
}
