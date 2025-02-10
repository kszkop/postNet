\name{ptn_model}
\alias{ptn_model}
\title{Retrieve the step-wise regression model generated during feature integration from the postNetData object}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptn_model(ptn, analysis_type, model, comparison)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{analysis_type}{
%%     ~~Describe \code{analysis_type} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{comparison}{
%%     ~~Describe \code{comparison} here~~
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
function (ptn, analysis_type, model, comparison) 
{
    check_ptn(ptn)
    check_analysis_type(analysis_type)
    check_model(model, analysis_type = analysis_type)
    if (!check_number(comparison)) {
        stop(paste("Please provide one numeric value for ", comparison, 
            sep = ""))
    }
    tmpIn <- slot(ptn@analysis@featureIntegration, analysis_type)
    if (comparison > length(tmpIn)) {
        stop(paste("There are only ", length(tmpIn), " comparisons", 
            sep = ""))
    }
    tmpIn <- tmpIn[[comparison]]
    tmpOut <- slot(tmpIn, model)
    return(tmpOut)
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
