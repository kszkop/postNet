\name{signaturesHeatmap}
\alias{signaturesHeatmap}
\title{Evaluate regulation of gene signatures in an anota2seq analysis with heatmaps
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
signaturesHeatmap(signatureList, ads, contrast, contrastNames, unit, RegMode, pdfName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signatureList}{
%%     ~~Describe \code{signatureList} here~~
}
  \item{ads}{
%%     ~~Describe \code{ads} here~~
}
  \item{contrast}{
%%     ~~Describe \code{contrast} here~~
}
  \item{contrastNames}{
%%     ~~Describe \code{contrastNames} here~~
}
  \item{unit}{
%%     ~~Describe \code{unit} here~~
}
  \item{RegMode}{
%%     ~~Describe \code{RegMode} here~~
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
function (signatureList, ads, contrast, contrastNames, unit, 
    RegMode, pdfName) 
{
    tableFinal <- matrix(NA, nrow = length(signatureList), ncol = length(contrasts))
    row.names(tableFinal) <- names(signatureList)
    colnames(tableFinal) <- contrastNames
    for (cont in 1:length(contrast)) {
        regData <- data.frame(geneSymb = rownames(ads@dataP))
        if (RegMode == "translation") {
            regData$ApvEff <- ads@translation@apvStatsRvm[[cont]][, 
                "apvEff"]
        }
        if (RegMode == "buffering") {
            regData$ApvEff <- ads@buffering@apvStatsRvm[[cont]][, 
                "apvEff"]
        }
        if (RegMode == "poly") {
            regData$ApvEff <- ads@translatedmRNA@apvStatsRvm[[cont]][, 
                "apvEff"]
        }
        if (RegMode == "total") {
            regData$ApvEff <- ads@totalmRNA@apvStatsRvm[[cont]][, 
                "apvEff"]
        }
        percOut <- as.numeric()
        if (unit == "FDR") {
            fdrOut <- as.numeric()
        }
        for (sign in 1:length(signatureList)) {
            regData[, 3] <- "bkg"
            regData[, 3][regData$geneSymb \%in\% signatureList[[sign]]] <- names(signatureList)[sign]
            colnames(regData)[3] <- "signature"
            tmpBg <- sort(as.numeric(regData[regData$signature == 
                "bkg", ]$ApvEff))
            ecdfBg <- 1:length(tmpBg)/length(tmpBg)
            if (unit == "FDR") {
                percentileBG <- tmpBg[which(ecdfBg >= 0.5)[1]]
            }
            else {
                percentileBG <- tmpBg[which(ecdfBg >= as.numeric(gsub("p", 
                  "0.", unit)))[1]]
            }
            tmpSign <- sort(as.numeric(regData[regData$signature == 
                names(signatureList)[sign], ]$ApvEff))
            ecdfSign <- 1:length(tmpSign)/length(tmpSign)
            if (unit == "FDR") {
                percentileSign <- tmpSign[which(ecdfSign >= 0.5)[1]]
            }
            else {
                percentileSign <- tmpSign[which(ecdfSign >= as.numeric(gsub("p", 
                  "0.", unit)))[1]]
            }
            percentileDiff <- percentileSign - percentileBG
            percOut[sign] <- percentileDiff
            if (unit == "FDR") {
                pval <- as.numeric(wilcox.test(as.numeric(regData[regData$signature == 
                  names(signatureList)[sign], ]$ApvEff), as.numeric(regData[regData$signature == 
                  "bkg", ]$ApvEff), alternative = "two.sided")[3])
                fdrOut[sign] <- pval
            }
        }
        if (unit == "FDR") {
            adjFDR <- -log10(p.adjust(fdrOut))
            fdrDirec <- ifelse(percOut >= 0, adjFDR * 1, adjFDR * 
                -1)
            tableFinal[, cont] <- fdrDirec
        }
        else {
            tableFinal[, cont] <- percOut
        }
    }
    if (unit == "FDR") {
        if (max(abs(c(min(as.vector(tableFinal)), max(as.vector(tableFinal))))) > 
            10) {
            breaks <- seq(-10, 10, length.out = 20)
            len <- length(breaks) - 1
        }
        else {
            breaks <- seq(min(as.vector(tableFinal)), max(as.vector(tableFinal)), 
                length.out = 25)
            len <- length(breaks) - 1
        }
        keyL <- paste("-log10", unit, sep = " ")
    }
    else {
        breaks <- seq(-max(abs(min(as.vector(tableFinal))), abs(max(as.vector(tableFinal)))), 
            max(abs(min(as.vector(tableFinal))), abs(max(as.vector(tableFinal)))), 
            length.out = 25)
        len <- length(breaks) - 1
        keyL <- unit
    }
    if (length(contrast) == 1) {
        tableFinal <- cbind(tableFinal, rep(0, nrow(tableFinal)))
    }
    pdf(ifelse(is.null(pdfName), "heatmap.pdf", paste(pdfName, 
        "heatmap.pdf", sep = "_")), width = 8, height = 8, useDingbats = F)
    par(mar = c(10, 5, 5, 10), bty = "l", font = 2, font.axis = 2, 
        font.lab = 2, cex.axis = 0.9, cex.main = 0.7, cex.lab = 0.9)
    gplots::heatmap.2(tableFinal, trace = "none", breaks = breaks, 
        Rowv = TRUE, col = rev(colorRampPalette(RColorBrewer::brewer.pal(name = "RdBu", 
            n = 11))(len)), key.xlab = keyL, key.title = "", 
        dendrogram = "row", density.info = "none", tracecol = NULL, 
        margins = c(10, 19), lhei = c(1, 6), lwid = c(0.5, 1), 
        cexRow = 0.9, cexCol = 0.9)
    dev.off()
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
