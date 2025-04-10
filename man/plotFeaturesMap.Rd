\name{plotFeaturesMap}
\alias{plotFeaturesMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotFeaturesMap(ptn, regOnly = TRUE, comparisons = NULL, featSel, scaled = TRUE, remExtreme = NULL, pdfName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{regOnly}{
%%     ~~Describe \code{regOnly} here~~
}
  \item{comparisons}{
%%     ~~Describe \code{comparisons} here~~
}
  \item{featSel}{
%%     ~~Describe \code{featSel} here~~
}
  \item{scaled}{
%%     ~~Describe \code{scaled} here~~
}
  \item{remExtreme}{
%%     ~~Describe \code{remExtreme} here~~
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
function (ptn, regOnly = TRUE, comparisons = NULL, featSel, scaled = TRUE, 
    remExtreme = NULL, pdfName = NULL) 
{
    check_ptn(ptn)
    if (!is.null(ptn_features(ptn))) {
        featuresTmp <- ptn_features(ptn)
    }
    else {
        stop("features analysis is null, please run featureIntegration analysis first")
    }
    if (!check_featSel(featSel, features = featuresTmp)) {
        stop("Please provide character vector of features to be plotted. The recommedation is to provide these that were selected in the final featureIntegration model.\n         These can be obtained by running ptn_selectedFeatures function and providing ptn object, analysis_type (ie lm or rf) and comparison , it would be a consecutivenumber from the comparisons provided running featureIntegration")
    }
    featuresTmp <- featuresTmp[, colnames(featuresTmp) \%in\% featSel]
    if (!check_logical(regOnly)) {
        stop("'regOnly' can only be TRUE or FALSE")
    }
    if (!check_logical(scaled)) {
        stop("'scaled' can only be TRUE or FALSE")
    }
    if (isTRUE(regOnly)) {
        if (!check_comparisons(comparisons)) {
            stop("'comparisons' must be a list of numeric vector for paired comparisons example: list(c(0,2),c(0,1)). 0 is always a background.")
        }
        if (length(which(unique(unlist(comparisons)) == 0)) > 
            0 && is.null(ptn_background(ptn))) {
            stop(" 0 is always a background, but no background provided")
        }
        if (length(comparisons) != 1) {
            stop("Although potentially posssible, please run each comparison separetely to fit well to your feature selection")
        }
    }
    if (!is.null(remExtreme)) {
        if (!check_number(remExtreme)) {
            stop("remExtreme should be a number indicating percentile of extreme values to be removed for colouring purposes ")
        }
        else {
            if (remExtreme <= 0 | remExtreme >= 1) {
                stop("remExtreme can be only 1 number between (0,1). It represents percentile in decimal of extreme values from both sides of distribution to be removed only for colouring scale")
            }
        }
    }
    if (length(which(unique(unlist(comparisons)) == 0)) > 0 && 
        is.null(ptn_background(ptn))) {
        stop(" 0 is always a background, but no background provided")
    }
    if (isTRUE(regOnly)) {
        resOut <- resQuant(qvec = ptn_effect(ptn), ptn = ptn)
        if (names(resOut)[1] == "background") {
            compTmp <- comparisons[[1]] + 1
        }
        else {
            compTmp <- comparisons[[1]]
        }
        listSel <- c(names(resOut[[compTmp[1]]]), names(resOut[[compTmp[2]]]))
        features <- featuresTmp[row.names(featuresTmp) \%in\% listSel, 
            ]
    }
    if (isTRUE(scaled)) {
        featuresOut <- scale(na.omit(features), center = T)
    }
    else {
        featuresOut <- na.omit(features)
    }
    fmapRes <- umap::umap(featuresOut, n_components = 2)
    fmapRes <- fmapResOut <- as.data.frame(fmapRes$layout)
    colnames(fmapRes) <- colnames(fmapResOut) <- c("fUMAP1", 
        "fUMAP2")
    fmapRes$Gene <- rownames(featuresOut)
    effTmp <- ptn_effect(ptn)
    eff <- effTmp[match(row.names(featuresOut), names(effTmp))]
    effect_fmap <- plot_fmap(fmapRes, colVec = eff, remExtreme = remExtreme, 
        name = "Effect")
    for (feat in featSel) {
        featTmp <- features[, feat]
        feature_fmap <- plot_fmap(fmapRes, colVec = featTmp, 
            remExtreme = remExtreme, name = feat)
        pdf(ifelse(is.null(pdfName), paste(feat, "_fumap.pdf", 
            sep = ""), paste(pdfName, feat, "fumap.pdf", sep = "_")), 
            width = 16, height = 8, useDingbats = F)
        par(mar = c(5, 5, 5, 5), bty = "l", font = 2, font.axis = 2, 
            font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.9)
        gridExtra::grid.arrange(effect_fmap$legend[[1]], effect_fmap$mainPlot, 
            feature_fmap$mainPlot, feature_fmap$legend[[1]], 
            ncol = 4, widths = c(1, 4, 4, 1))
        dev.off()
    }
    ptn@analysis@featureIntegration@featureMap <- fmapResOut
    return(ptn)
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
