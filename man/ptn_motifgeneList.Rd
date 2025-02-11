\name{ptn_motifgeneList}
\alias{ptn_motifgeneList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ptn_motifgeneList
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptn_motifgeneList(x, region, geneList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{region}{
%%     ~~Describe \code{region} here~~
}
  \item{geneList}{
%%     ~~Describe \code{geneList} here~~
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
new("standardGeneric", .Data = function (x, region, geneList) 
standardGeneric("ptn_motifgeneList"), generic = structure("ptn_motifgeneList", package = "postNet"), 
    package = "postNet", group = list(), valueClass = character(0), 
    signature = c("x", "region", "geneList"), default = NULL, 
    skeleton = (function (x, region, geneList) 
    stop(gettextf("invalid call in method dispatch to '\%s' (no default method)", 
        "ptn_motifgeneList"), domain = NA))(x, region, geneList))
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
