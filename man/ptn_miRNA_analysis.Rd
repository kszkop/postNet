\name{ptn_miRNA_analysis}
\alias{ptn_miRNA_analysis}
\title{Retrieve the results of miRNA analysis from the postNetData object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptn_miRNA_analysis(ptn, direction, threshold)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
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
function (ptn, direction, threshold) 
{
    check_direction(tolower(direction))
    check_ptn(ptn)
    if (!check_number(threshold)) {
        stop(paste("Please provide one numeric value for ", threshold, 
            sep = ""))
    }
    if (is.null(slot(ptn@analysis, "miRNA"))) {
        stop("Please run miRNAanalysis first")
    }
    else {
        miRNAres <- ptn@analysis@miRNA@miRNA_analysis
    }
    if (tolower(direction) == "greater") {
        resOut <- miRNAres$greater
    }
    else if (tolower(direction) == "less") {
        resOut <- miRNAres$less
    }
    resOut <- resOut[which(resOut[, 4] < threshold), ]
    if (nrow(resOut) > 0) {
        resOut <- resOut[, c(1, 2, 5, 3, 4)]
        resOut <- data.frame(id = row.names(resOut), resOut, 
            row.names = NULL)
    }
    else {
        message("there are no miRNAs to output")
    }
    return(resOut)
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
