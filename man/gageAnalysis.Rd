\name{gageAnalysis}
\alias{gageAnalysis}
\title{Perform GAGE analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gageAnalysis(ptn, category, genesSlopeFiltOut = NULL, maxSize = 500, minSize = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{category}{
%%     ~~Describe \code{category} here~~
}
  \item{genesSlopeFiltOut}{
%%     ~~Describe \code{genesSlopeFiltOut} here~~
}
  \item{maxSize}{
%%     ~~Describe \code{maxSize} here~~
}
  \item{minSize}{
%%     ~~Describe \code{minSize} here~~
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
function (ptn, category, genesSlopeFiltOut = NULL, maxSize = 500, 
    minSize = 10) 
{
    check_ptn(ptn)
    species <- ptn_species(ptn)
    if (!species \%in\% c("human", "mouse")) {
        stop("This option is only  available for species: human and mouse at the moment")
    }
    if (!check_number(maxSize) | !check_number(minSize)) {
        stop("please provide numeric value")
    }
    if (minSize <= 0 | maxSize <= 0) {
        stop("size parameters must be positive")
    }
    if (maxSize <= minSize) {
        stop("maxSize must be greater than minSize")
    }
    gageOut <- list()
    effTmp <- ptn_effect(ptn)
    if (!is.null(genesSlopeFiltOut)) {
        effIn <- effTmp[!names(effTmp) \%in\% genesSlopeFiltOut]
    }
    else {
        effIn <- effTmp
    }
    convTab <- convertSymbolToEntrezID(names(effIn), species = species)
    names(convTab) <- convertEntrezIDToSymbol(convTab, species = species)
    names(effIn) <- as.character(convTab)[match(names(effIn), 
        names(convTab))]
    effIn <- effIn[!is.na(names(effIn))]
    rankIn <- effIn[order(effIn, decreasing = T)]
    GAGEout <- new("postNetGAGE", BP = NULL, CC = NULL, MF = NULL, 
        KEGG = NULL)
    for (sel in category) {
        resOut <- list()
        if (sel == "KEGG") {
            pathwayKegg <- gage::kegg.gsets(species = species, 
                id.type = "entrez")
            pathwaysIn <- pathwayKegg$kg.sets
        }
        else {
            if (!exists("pathwayGO")) {
                pathwayGO <- gage::go.gsets(species = species, 
                  id.type = "EG")
            }
            if (sel == "BP") {
                pathwaysIn <- pathwayGO$go.sets[pathwayGO$go.subs$BP]
            }
            else if (sel == "MF") {
                pathwaysIn <- pathwayGO$go.sets[pathwayGO$go.subs$MF]
            }
            else if (sel == "CC") {
                pathwaysIn <- pathwayGO$go.sets[pathwayGO$go.subs$CC]
            }
        }
        resOut <- gage::gage(exprs = rankIn, gsets = pathwaysIn, 
            set.size = c(minSize, maxSize), rank.test = TRUE, 
            use.fold = FALSE)
        pathwaysGenes <- list()
        for (pw in 1:length(pathwaysIn)) {
            glist <- unlist(pathwaysIn[pw])
            tmpGenes <- glist[glist \%in\% names(rankIn)]
            if (length(tmpGenes) > 0) {
                tmpGenes_conv <- names(convTab)[match(as.character(tmpGenes), 
                  convTab)]
            }
            else {
                tmpGenes_conv <- NA
            }
            pathwaysGenes[[pw]] <- paste(tmpGenes_conv, collapse = ":")
        }
        names(pathwaysGenes) <- names(pathwaysIn)
        pathwaysGenes <- pathwaysGenes[as.numeric(which(unlist(lapply(pathwaysGenes, 
            function(x) x != "NA"))))]
        if (nrow(resOut$greater) > 0) {
            grTmpG <- resOut$greater
            grTmpG <- grTmpG[, c(1, 2, 5, 3, 4)]
            grTmpG <- data.frame(id = row.names(grTmpG), grTmpG, 
                row.names = NULL)
            grTmpG$Genes <- as.character(pathwaysGenes)[match(grTmpG$id, 
                names(pathwaysGenes))]
            resOut$greater <- grTmpG
        }
        if (nrow(resOut$less) > 0) {
            grTmpL <- resOut$less
            grTmpL <- grTmpL[, c(1, 2, 5, 3, 4)]
            grTmpL <- data.frame(id = row.names(grTmpL), grTmpL, 
                row.names = NULL)
            grTmpL$Genes <- as.character(pathwaysGenes)[match(grTmpL$id, 
                names(pathwaysGenes))]
            resOut$less <- grTmpL
        }
        slot(GAGEout, sel) <- resOut
    }
    ptn@analysis@GAGE <- GAGEout
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
