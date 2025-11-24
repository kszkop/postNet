\name{get_signatures}
\alias{get_signatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_signatures(species)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species}{
%%     ~~Describe \code{species} here~~
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
function (species) 
{
    if (!is_valid_species(species)) {
        stop("Please specify a species, at the moment only 'human' or 'mouse' are available).")
    }
    currTmp <- list.files(system.file("extdata/signatures", package = "postNetParcel"))
    if (!species \%in\% currTmp) {
        stop("This option is currently only available for species 'human' and 'mouse'. Please use the options 'custom' and 'customFile' to provide annotations for other species.")
    }
    signatures <- readRDS(system.file(paste("extdata/signatures", 
        species, sep = "/"), paste(species, "Signatures.rds", 
        sep = ""), package = "postNetParcel"))
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
