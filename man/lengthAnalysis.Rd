\name{lengthAnalysis}
\alias{lengthAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparisons of lengths
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lengthAnalysis(ptn, region, comparisons = NULL, plotOut = TRUE, plotType = "boxplot", pdfName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{region}{
%%     ~~Describe \code{region} here~~
}
  \item{comparisons}{
%%     ~~Describe \code{comparisons} here~~
}
  \item{plotOut}{
%%     ~~Describe \code{plotOut} here~~
}
  \item{plotType}{
%%     ~~Describe \code{plotType} here~~
}
  \item{pdfName}{
%%     ~~Describe \code{pdfName} here~~
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
function (ptn, region, comparisons = NULL, plotOut = TRUE, plotType = "boxplot", 
    pdfName = NULL) 
{
    check_ptn(ptn)
    check_region(region)
    if (!check_logical(plotOut)) {
        stop("'plotOut' can only be only be logical: TRUE of FALSE ")
    }
    if (isTRUE(plotOut)) {
        if (!is.null(plotType)) {
            check_plotType(plotType)
        }
        else {
            stop("Please provide 'plotType' to select option for plotting, from: 'boxplot','violin ,'ecdf'. ")
        }
    }
    if (!is.null(comparisons)) {
        if (!check_comparisons(comparisons)) {
            stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
        }
        if (length(which(unique(unlist(comparisons)) == 0)) > 
            0 && is.null(ptn_background(ptn))) {
            stop(" 0 is always a background, but no background provided")
        }
    }
    lengthFinal <- list()
    for (reg in region) {
        seqTmp <- ptn_sequences(ptn, region = reg)
        lenTmp <- as.numeric(sapply(seqTmp, function(x) length(seqinr::s2c(x))))
        lenForAnalysis <- log2(as.numeric(lenTmp))
        names(lenForAnalysis) <- ptn_geneID(ptn, region = reg)
        if (isTRUE(plotOut)) {
            resOut <- resQuant(qvec = lenForAnalysis, ptn = ptn)
            if (length(resOut) == 0) {
                stop("There are no regulated genes. Check the input or run without indicating regulation and comparisons")
            }
            colOut <- colPlot(ptn)
            pdf(ifelse(is.null(pdfName), paste(reg, plotType, 
                "lengthAnalysis.pdf", sep = "_"), paste(pdfName, 
                reg, plotType, "lengthAnalysis.pdf", sep = "_")), 
                width = 8, height = 8, useDingbats = F)
            ylabel = "Length (Log2 scale)"
            plotUtils(resOut, colOut, comparisons, ylabel = ylabel, 
                plotType = plotType)
            dev.off()
        }
        lengthFinal[[paste(reg, "length", sep = "_")]] <- lenForAnalysis
    }
    return(lengthFinal)
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
