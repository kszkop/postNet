\name{ptn_GO}
\alias{ptn_GO}
\title{Retrieve results of GO term analysis from the postNetData object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptn_GO(ptn, category, geneList, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{category}{
%%     ~~Describe \code{category} here~~
}
  \item{geneList}{
%%     ~~Describe \code{geneList} here~~
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
function (ptn, category, geneList, threshold) 
{
    check_ptn(ptn)
    check_category(category)
    if (length(category) != 1) {
        stop("Please provide only one category")
    }
    if (!check_number(threshold)) {
        stop(paste("Please provide one numeric value for ", threshold, 
            sep = ""))
    }
    if (!any(geneList \%in\% names(ptn_geneList(ptn)))) {
        stop("None of the regulatory geneList in ptn")
    }
    if (is.null(slot(ptn@analysis, "GO"))) {
        stop("Please run GO analysis first")
    }
    else {
        GOres <- slot(ptn@analysis@GO, category)
    }
    GOresOut <- GOres[[which(geneList == names(GOres))]]@result
    GOresOut <- GOresOut[which(GOresOut[, 6] < threshold), ]
    if (nrow(GOresOut) > 0) {
        GOresOut <- data.frame(GOresOut, row.names = NULL)
    }
    else {
        message("there are no GO categories to output")
    }
    return(GOresOut)
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
