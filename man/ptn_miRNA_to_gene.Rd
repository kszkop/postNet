\name{ptn_miRNA_to_gene}
\alias{ptn_miRNA_to_gene}
\title{Retrieve other gene targets of miRNA from the postNetData object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptn_miRNA_to_gene(ptn, miRNAs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{miRNAs}{
%%     ~~Describe \code{miRNAs} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ptn, miRNAs) 
{
    check_ptn(ptn)
    if (is.null(slot(ptn@analysis, "miRNA_analysis"))) {
        stop("Please run miRNAanalysis first")
    }
    else {
        miRNATmp <- ptn@analysis@miRNA@miRNA_to_gene
    }
    miRNAsOut <- miRNATmp[which(names(miRNATmp) \%in\% miRNAs)]
    return(miRNAsOut)
  }
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
