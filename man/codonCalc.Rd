\name{codonCalc}
\alias{codonCalc}
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
codonCalc(ptn, featsel, analysis = "codon", unit = "count", comparisons = NULL, plotOut = TRUE, plotType = "ecdf", pdfName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{featsel}{
%%     ~~Describe \code{featsel} here~~
}
  \item{analysis}{
%%     ~~Describe \code{analysis} here~~
}
  \item{unit}{
%%     ~~Describe \code{unit} here~~
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
##--	or standard data sets, see data().

## The function is currently defined as
function (ptn, featsel, analysis = "codon", unit = "count", comparisons = NULL, 
    plotOut = TRUE, plotType = "ecdf", pdfName = NULL) 
{
    check_ptn(ptn)
    if (!is.null(ptn_codonAnalysis(ptn))) {
        codonsAll <- ptn_codonAnalysis(ptn)
        check_codonIn(codonsAll)
    }
    else {
        stop("codons analysis is null, please provide valid 'postNetData' object and run codonUsage analysis")
    }
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
    if (!is_valid_analysis(analysis)) {
        stop("Please provide an 'analysis'. It can only be 'codon' or 'AA'")
    }
    if (analysis == "codon") {
        if (!check_codons(featsel)) {
            stop("featsel does not contain valid codons")
        }
    }
    if (!isUnit(unit)) {
        stop("'unit' must be one from these: 'count' or 'freq'")
    }
    codonCalcOut <- list()
    for (i in 1:length(featsel)) {
        featTmp <- featsel[[i]]
        featNameTmp <- names(featsel)[i]
        if (featNameTmp == 0) {
            stop("'featSel should be a named list")
        }
        nameTmp <- ifelse(is.null(pdfName), paste("features", 
            i, "codonCalc.pdf", sep = "_"), paste(pdfName, "features", 
            i, "codonCalc.pdf", sep = "_"))
        nameOut <- nameTmp
        if (analysis == "codon") {
            codonTmp <- codonsAll[codonsAll$codon \%in\% featTmp, 
                ]
        }
        else if (analysis == "AA") {
            codonTmp <- codonsAll[codonsAll$AA \%in\% featTmp, 
                ]
        }
        if (unit == "count") {
            tmp <- codonTmp \%>\% group_by(geneID) \%>\% summarise(count = sum(count))
            codonCalcOutTmp <- tmp$count
        }
        else if (unit == "freq") {
            tmp <- codonTmp \%>\% group_by(geneID) \%>\% summarise(freq = sum(frequency))
            codonCalcOutTmp <- tmp$freq
        }
        names(codonCalcOutTmp) <- tmp$geneID
        codonCalcOut[[featNameTmp]] <- codonCalcOutTmp
        if (isTRUE(plotOut)) {
            resOut <- resQuant(qvec = codonCalcOutTmp, ptn = ptn)
            if (length(resOut) == 0) {
                stop("There are no regulated genes. Check the input or run without indicating regulation and comparisons")
            }
            colOut <- colPlot(ptn)
            pdf(nameOut, width = 8, height = 8, useDingbats = F)
            ylabel <- paste("codon usage(", unit, ")", sep = "")
            plotPostNet(resOut, colOut, comparisons, ylabel = ylabel, 
                plotType = plotType)
            dev.off()
        }
    }
    return(codonCalcOut)
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
