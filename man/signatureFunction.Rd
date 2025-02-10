\name{signatureFunction}
\alias{signatureFunction}
\title{Evaluate regulation of gene signatures in an anota2seq analysis 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
signatureFunction(signatureList, generalName, dataName, effects_names = c("total mRNA log2FC", "polysome associated mRNA log2FC", "buffering log2FC", "translation log2FC"), colours, ads, contrast, xlim = NULL, scatterXY = NULL, tableCex, pdfName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signatureList}{
%%     ~~Describe \code{signatureList} here~~
}
  \item{generalName}{
%%     ~~Describe \code{generalName} here~~
}
  \item{dataName}{
%%     ~~Describe \code{dataName} here~~
}
  \item{effects_names}{
%%     ~~Describe \code{effects_names} here~~
}
  \item{colours}{
%%     ~~Describe \code{colours} here~~
}
  \item{ads}{
%%     ~~Describe \code{ads} here~~
}
  \item{contrast}{
%%     ~~Describe \code{contrast} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{scatterXY}{
%%     ~~Describe \code{scatterXY} here~~
}
  \item{tableCex}{
%%     ~~Describe \code{tableCex} here~~
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
function (signatureList, generalName, dataName, effects_names = c("total mRNA log2FC", 
    "polysome associated mRNA log2FC", "buffering log2FC", "translation log2FC"), 
    colours, ads, contrast, xlim = NULL, scatterXY = NULL, tableCex, 
    pdfName) 
{
    signNames <- names(signatureList)
    regData <- data.frame(geneSymb = rownames(ads@dataP))
    regData$totalApvEff <- ads@totalmRNA@apvStatsRvm[[contrast]][, 
        "apvEff"]
    regData$polyApvEff <- ads@translatedmRNA@apvStatsRvm[[contrast]][, 
        "apvEff"]
    regData$buffApvEff <- ads@buffering@apvStatsRvm[[contrast]][, 
        "apvEff"]
    regData$translationApvEff <- ads@translation@apvStatsRvm[[contrast]][, 
        "apvEff"]
    if (is.null(scatterXY)) {
        scatterXY <- roundNice(max(abs(c(range(regData$totalApvEff), 
            range(regData$polyApvEff)))), direction = "up")
    }
    pdf(ifelse(is.null(pdfName), paste("data_", dataName, "_signature_", 
        generalName, ".pdf", sep = ""), paste(pdfName, paste("data_", 
        dataName, "_signature_", generalName, ".pdf", sep = ""), 
        sep = "_")), width = 18, height = 4, useDingbats = F)
    par(mfrow = c(1, 5), mar = c(5, 5, 6, 4), bty = "l", font = 2, 
        font.axis = 2, font.lab = 2, cex.axis = 1.25, cex.main = 1.9, 
        cex.lab = 2)
    plot(regData$totalApvEff, regData$polyApvEff, pch = 16, cex = 1.9, 
        col = "grey75", ylab = effects_names[2], xlab = effects_names[1], 
        main = paste(paste("Data: ", dataName, sep = ""), "\n", 
            paste("Signature: ", generalName, ""), sep = ""), 
        xlim = c(-scatterXY, scatterXY), ylim = c(-scatterXY, 
            scatterXY))
    abline(v = 0)
    abline(h = 0)
    if (any(duplicated(unlist(signatureList)))) {
        cat("There are some genes that overlap between signatures. Separate background for each will be used")
        for (i in 1:length(signatureList)) {
            regData[, 5 + i] <- "bkg"
            regData[, 5 + i][regData$geneSymb \%in\% signatureList[[i]]] <- names(signatureList)[i]
            points(regData$totalApvEff[regData[, 5 + i] == names(signatureList)[i]], 
                regData$polyApvEff[regData[, 5 + i] == names(signatureList)[i]], 
                col = colours[i], pch = 16, cex = 1.7)
        }
    }
    else {
        regData$signature <- "bkg"
        for (i in 1:length(signatureList)) {
            regData$signature[regData$geneSymb \%in\% signatureList[[i]]] <- names(signatureList)[i]
            points(regData$totalApvEff[regData$signature == names(signatureList)[i]], 
                regData$polyApvEff[regData$signature == names(signatureList)[i]], 
                col = colours[i], pch = 16, cex = 1.7)
        }
    }
    for (eff in 5:2) {
        tableOut <- matrix(NA, nrow = length(signNames), ncol = 5)
        colnames(tableOut) <- c("signature", "Wilcox_pval", "q25", 
            "q50", "q75")
        tableOut[, 1] <- signNames
        if (any(duplicated(unlist(signatureList)))) {
            tmpBgOut <- list()
            for (i in 1:length(signatureList)) {
                tmpBg <- sort(as.numeric(regData[regData[, (5 + 
                  i)] == "bkg", ][, eff]))
                ecdfBg <- 1:length(tmpBg)/length(tmpBg)
                bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
                bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
                bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
                tmpBgOut[[i]] <- tmpBg
                tableOut[i, 2] <- format(as.numeric(wilcox.test(as.numeric(regData[regData[, 
                  (5 + i)] == signNames[i], ][, eff]), as.numeric(regData[regData[, 
                  (5 + i)] == "bkg", ][, eff]), alternative = "two.sided")[3]), 
                  scientific = TRUE, digits = 2)
                tmpSign <- sort(as.numeric(regData[regData[, 
                  (5 + i)] == signNames[i], ][, eff]))
                ecdfSign <- 1:length(tmpSign)/length(tmpSign)
                tableOut[i, 3] <- format(tmpSign[which(ecdfSign >= 
                  0.25)[1]] - bg_025, digits = 2)
                tableOut[i, 4] <- format(tmpSign[which(ecdfSign >= 
                  0.5)[1]] - bg_05, digits = 2)
                tableOut[i, 5] <- format(tmpSign[which(ecdfSign >= 
                  0.75)[1]] - bg_075, digits = 2)
            }
            if (!is.null(xlim)) {
                xmin <- xlim[1]
                xmax <- xlim[2]
            }
            else {
                xmin <- ifelse(as.numeric(quantile(as.numeric(unlist(tmpBgOut)), 
                  0.01)) < 0, -roundNice(abs(as.numeric(quantile(as.numeric(tmpBg), 
                  0.01))), direction = "up"), roundNice(as.numeric(quantile(as.numeric(unlist(tmpBgOut)), 
                  0.01))), direction = "up")
                xmax <- ifelse(as.numeric(quantile(as.numeric(unlist(tmpBgOut)), 
                  0.99)) < 0, -roundNice(abs(as.numeric(quantile(as.numeric(tmpBg), 
                  0.01))), direction = "up"), roundNice(as.numeric(quantile(as.numeric(unlist(tmpBgOut)), 
                  0.99))), direction = "up")
            }
            plot(ecdf(as.numeric(regData[regData[, (5 + i)] == 
                "bkg", ][, eff])), col = "grey55", main = "", 
                xlab = effects_names[eff - 1], verticals = TRUE, 
                do.p = FALSE, lwd = 3, xlim = c(xmin, xmax))
            legend(xmin, 0.95, fill = "grey55", border = "grey55", 
                "Background", bty = "n", cex = 1.3)
            for (i in 1:length(signatureList)) {
                lines(ecdf(as.numeric(regData[regData[, (5 + 
                  i)] == signNames[i], ][, eff])), col = colours[i], 
                  main = "", xlab = "", verticals = TRUE, do.p = FALSE, 
                  lwd = 3)
            }
            plotrix::addtable2plot(xmin - abs((xmin * 0.1)), 
                1.05, tableOut, bty = "n", display.rownames = FALSE, 
                hlines = FALSE, vlines = TRUE, title = "", cex = tableCex, 
                bg = colours, xpad = 0.2, ypad = 1.4)
        }
        else {
            tmpBg <- sort(as.numeric(regData[regData$signature == 
                "bkg", ][, eff]))
            ecdfBg <- 1:length(tmpBg)/length(tmpBg)
            bg_025 <- tmpBg[which(ecdfBg >= 0.25)[1]]
            bg_05 <- tmpBg[which(ecdfBg >= 0.5)[1]]
            bg_075 <- tmpBg[which(ecdfBg >= 0.75)[1]]
            for (i in 1:length(signatureList)) {
                tableOut[i, 2] <- format(as.numeric(wilcox.test(as.numeric(regData[regData$signature == 
                  signNames[i], ][, eff]), as.numeric(regData[regData$signature == 
                  "bkg", ][, eff]), alternative = "two.sided")[3]), 
                  scientific = TRUE, digits = 2)
                tmpSign <- sort(as.numeric(regData[regData$signature == 
                  signNames[i], ][, eff]))
                ecdfSign <- 1:length(tmpSign)/length(tmpSign)
                tableOut[i, 3] <- format(tmpSign[which(ecdfSign >= 
                  0.25)[1]] - bg_025, digits = 2)
                tableOut[i, 4] <- format(tmpSign[which(ecdfSign >= 
                  0.5)[1]] - bg_05, digits = 2)
                tableOut[i, 5] <- format(tmpSign[which(ecdfSign >= 
                  0.75)[1]] - bg_075, digits = 2)
            }
            if (!is.null(xlim)) {
                xmin <- xlim[1]
                xmax <- xlim[2]
            }
            else {
                xmin <- ifelse(as.numeric(quantile(as.numeric(tmpBg), 
                  0.01)) < 0, -roundNice(abs(as.numeric(quantile(as.numeric(tmpBg), 
                  0.01))), direction = "up"), roundNice(as.numeric(quantile(as.numeric(tmpBg), 
                  0.01))), direction = "up")
                xmax <- ifelse(as.numeric(quantile(as.numeric(tmpBg), 
                  0.99)) < 0, -roundNice(abs(as.numeric(quantile(as.numeric(tmpBg), 
                  0.01))), direction = "up"), roundNice(as.numeric(quantile(as.numeric(tmpBg), 
                  0.99))), direction = "up")
            }
            plot(ecdf(as.numeric(regData[regData$signature == 
                "bkg", ][, eff])), col = "grey55", main = "", 
                xlab = effects_names[eff - 1], verticals = TRUE, 
                do.p = FALSE, lwd = 3, xlim = c(xmin, xmax))
            legend(xmin, 0.95, fill = "grey55", border = "grey55", 
                "Background", bty = "n", cex = 1.3)
            for (i in 1:length(signatureList)) {
                lines(ecdf(as.numeric(regData[regData$signature == 
                  signNames[i], ][, eff])), col = colours[i], 
                  main = "", xlab = "", verticals = TRUE, do.p = FALSE, 
                  lwd = 3)
            }
            plotrix::addtable2plot(xmin - abs((xmin * 0.1)), 
                1.05, tableOut, bty = "n", display.rownames = FALSE, 
                hlines = FALSE, vlines = TRUE, title = "", cex = tableCex, 
                bg = colours, xpad = 0.2, ypad = 1.4)
        }
    }
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
