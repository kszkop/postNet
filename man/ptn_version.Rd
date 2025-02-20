\docType{methods}
\name{ptn_version}
\alias{ptn_version,postNetData-method}
\title{Accessor for the RefSeq annotation release version from the "version" slot of a \code{postNetData} object
}

\description{The "version" slot holds the reslease version number of the RefSeq annotations used to generate the reference
sequences. Currently, this only applies when \code{source = 'load'} in the \link[postNet]{postNetStart} function. 
}

\usage{
\S4method{ptn_version}{postNetData}(x)
}

\arguments{
  \item{x}{a \code{postNetData} object.
  }
}

\value{A character specifying the release version of the RefSeq annotations used to generate the reference sequences.
}

\examples{

# load and create example data:
 data(postNetExample)

 myGenes <- postNetExample$geneList
 myBg <- postNetExample$background
 myEffect <- postNetExample$effect

# Initialize the postNetData object with custom gene lists:
 ptn <- postNetStart(
   geneList = myGenes,
   geneListcolours = c("#7FB7BE","#DB7F67"),
   customBg = myBg,
   effectMeasure = myEffect,
   selection = "random",
   source = "load", 
   species = "human"
 )

# Access version of the RefSeq annotation used to create the reference sequences
  ptn_version(ptn)

}

\seealso{
\code{\link[postNet]{postNetStart}} for details on selecting the appropriate reference sequence annotation versions.
}

