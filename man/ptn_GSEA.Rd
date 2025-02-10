\name{ptn_GSEA}
\alias{ptn_GSEA}
\title{Retrieve the results of GSEA analysis from the postNetData object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptn_GSEA(ptn, threshold = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
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
function (ptn, threshold = NULL) 
{
    check_ptn(ptn)
    if (is.null(slot(ptn@analysis, "GSEA"))) {
        stop("Please run GSEA analysis first")
    }
    else {
        gseaOut <- slot(ptn@analysis, "GSEA")
    }
    if (!is.null(threshold)) {
        if (!check_number(threshold)) {
            stop(paste("Please provide one numeric value for ", 
                threshold, sep = ""))
        }
        gseaOut <- gseaOut[which(gseaOut[, 8] < threshold), ]
    }
    if (nrow(gseaOut) == 0) {
        message("there are no GSEA gene sets to output")
    }
    return(gseaOut)
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
