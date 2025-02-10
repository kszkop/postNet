\name{ptn_GAGE}
\alias{ptn_GAGE}
\title{Retrieve the results of GAGE analysis from the postNetData object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptn_GAGE(ptn, category, direction, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{category}{
%%     ~~Describe \code{category} here~~
}
  \item{direction}{
%%     ~~Describe \code{direction} here~~
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
function (ptn, category, direction, threshold) 
{
    check_direction(tolower(direction))
    check_ptn(ptn)
    check_category(category)
    if (!check_number(threshold)) {
        stop(paste("Please provide one numeric value for ", threshold, 
            sep = ""))
    }
    if (length(category) != 1) {
        stop("Please provide only one category")
    }
    if (is.null(slot(ptn@analysis, "GAGE"))) {
        stop("Please run GAGE first ")
    }
    else {
        GAGEres <- slot(ptn@analysis@GAGE, category)
    }
    if (tolower(direction) == "greater") {
        resOut <- GAGEres$greater
    }
    else if (tolower(direction) == "less") {
        resOut <- GAGEres$less
    }
    resOut <- resOut[which(resOut[, 6] < threshold), ]
    if (nrow(resOut) == 0) {
        message("there are no terms enriched")
    }
    else {
        return(resOut)
    }
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
