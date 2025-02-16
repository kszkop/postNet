\name{ptn_sequences}
\alias{ptn_sequences}
\title{Retreive the reference annotation sequences from a postNetData obejct
}
\description{S4 method to retrieve the reference annotation sequences from a \code{postNetData} object.
}
\usage{
ptn_sequences(x, region)
}

\arguments{
  \item{x}{S4 object of class \code{postNetData}.
}
  \item{region}{The sequence region to be accessed. Can be either: \code{"UTR5"}, \code{"CDS"}, or
  \code{"UTR3"}.
  }
}

\value{
A character vector containing the reference annotation sequences corresponding to the specified region.
}

\seealso{

}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
new("standardGeneric", .Data = function (x, region) 
standardGeneric("ptn_sequences"), generic = structure("ptn_sequences", package = "postNet"), 
    package = "postNet", group = list(), valueClass = character(0), 
    signature = c("x", "region"), default = NULL, skeleton = (function (x, 
        region) 
    stop(gettextf("invalid call in method dispatch to '\%s' (no default method)", 
        "ptn_sequences"), domain = NA))(x, region))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
