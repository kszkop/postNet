\name{ptn_check_comparisons}
\alias{ptn_check_comparisons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptn_check_comparisons(ptn, analysis_type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{analysis_type}{
%%     ~~Describe \code{analysis_type} here~~
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
##--	or standard data sets, see data().

## The function is currently defined as
function (ptn, analysis_type) 
{
    check_ptn(ptn)
    check_analysis_type(analysis_type)
    tmpIn <- slot(ptn@analysis@featureIntegration, analysis_type)
    if (is.null(tmpIn)) {
        stop(paste("Please run ", analysis_type, "analysis first", 
            sep = ""))
    }
    else {
        print(names(tmpIn))
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
