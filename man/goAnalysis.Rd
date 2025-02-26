\name{goAnalysis}
\alias{goAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
goAnalysis(ptn, genesSlopeFiltOut = NULL, category, maxSize = 500, minSize = 10, counts = 10, FDR = 0.15, name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptn}{
%%     ~~Describe \code{ptn} here~~
}
  \item{genesSlopeFiltOut}{
%%     ~~Describe \code{genesSlopeFiltOut} here~~
}
  \item{category}{
%%     ~~Describe \code{category} here~~
}
  \item{maxSize}{
%%     ~~Describe \code{maxSize} here~~
}
  \item{minSize}{
%%     ~~Describe \code{minSize} here~~
}
  \item{counts}{
%%     ~~Describe \code{counts} here~~
}
  \item{FDR}{
%%     ~~Describe \code{FDR} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
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
function (ptn, genesSlopeFiltOut = NULL, category, maxSize = 500, 
    minSize = 10, counts = 10, FDR = 0.15, name = NULL) 
{
    check_ptn(ptn)
    species <- ptn_species(ptn)
    if (!species \%in\% c("human", "mouse")) {
        stop("This option is only  available for species: human and mouse at the moment")
    }
    check_category(category)
    if (!check_number(maxSize) | !check_number(minSize) | !check_number(counts) | 
        !check_number(FDR)) {
        stop("please provide numeric value")
    }
    if (minSize <= 0 | maxSize <= 0 | counts <= 0 | FDR < 0) {
        stop("size, counts and FDR parameters must be positive")
    }
    if (maxSize <= minSize) {
        stop("maxSize must be greater than minSize")
    }
    GOout <- list()
    res <- ptn_geneList(ptn)
    bg <- unlist(ptn_background(ptn))
    if (length(setdiff(bg, unique(unlist(res)))) == 0) {
        warning("Background seems not right as all genes are regulated")
    }
    if (!is.null(genesSlopeFiltOut)) {
        bg <- bg[!bg \%in\% genesSlopeFiltOut]
        res <- lapply(res, function(x) x[!x \%in\% genesSlopeFiltOut])
    }
    bg_entrezID <- convertSymbolToEntrezID(geneList = bg, species = species)
    res_entrezID <- lapply(res, function(x) convertSymbolToEntrezID(geneList = x, 
        species = species))
    GoLists <- res_entrezID[!sapply(res_entrezID, is.null)]
    GOout <- new("postNetGO", BP = NULL, CC = NULL, MF = NULL, 
        KEGG = NULL)
    for (sel in category) {
        resOut <- list()
        if (sel == "KEGG") {
            for (i in 1:length(GoLists)) {
                resTmp <- clusterProfiler::enrichKEGG(gene = GoLists[[i]], 
                  universe = bg_entrezID, organism = ifelse(species == 
                    "human", "hsa", ifelse(species == "mouse", 
                    "mmu", "ups")), pAdjustMethod = "BH", pvalueCutoff = 1, 
                  qvalueCutoff = 1, minGSSize = minSize, maxGSSize = maxSize)
                resOut[[names(GoLists)[i]]] <- resTmp
            }
        }
        else if (sel == "BP" | sel == "MF" | sel == "CC") {
            for (i in 1:length(GoLists)) {
                resTmp <- clusterProfiler::enrichGO(gene = GoLists[[i]], 
                  universe = bg_entrezID, OrgDb = ifelse(species == 
                    "human", "org.Hs.eg.db", ifelse(species == 
                    "mouse", "org.Mm.eg.db", "ups")), ont = sel, 
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, 
                  minGSSize = minSize, maxGSSize = maxSize, readable = TRUE, 
                  pool = FALSE)
                resOut[[names(GoLists)[i]]] <- resTmp
            }
        }
        for (i in 1:length(resOut)) {
            tabTmp <- resOut[[i]]@result
            tabTmp <- tabTmp[tabTmp$Count > counts, ]
            if (nrow(tabTmp) > 0) {
                tabTmp$p.adjust <- stats::p.adjust(tabTmp$pvalue, 
                  method = "BH")
            }
            tabTmp <- tabTmp[tabTmp$p.adjust < FDR, ]
            if (nrow(tabTmp) > 0) {
                geneIDs_temp <- tabTmp$geneID
                checkID <- check_id_type(seqinr::c2s(strsplit(geneIDs_temp[[1]], 
                  "/")[[1]][1:5]))
                if (checkID == "entrezID") {
                  tabTmp$geneID <- sapply(geneIDs_temp, function(x) paste(sort(convertEntrezIDToSymbol(unlist(strsplit(x, 
                    "/")), species = species)), collapse = ":"), 
                    USE.NAMES = F)
                }
                else {
                  tabTmp$geneID <- sapply(geneIDs_temp, function(x) paste(sort(unlist(strsplit(x, 
                    "/"))), collapse = ":"), USE.NAMES = F)
                }
            }
            else {
                message("No significant results")
            }
            tabTmp$Size <- as.numeric(sub("\\/.*", "", tabTmp$BgRatio))
            tabTmp <- tabTmp[, c(1, 2, 9, 10, 5, 6, 8)]
            resOut[[i]]@result <- tabTmp
        }
        resWrite <- lapply(resOut, function(x) x@result)
        nameOut <- ifelse(is.null(name), paste("GO_", sel, ".xlsx", 
            sep = ""), paste(name, "_GO_", sel, ".xlsx", sep = ""))
        WriteXLS::WriteXLS(resWrite, SheetNames = names(resWrite), 
            ExcelFileName = nameOut, row.names = FALSE)
        slot(GOout, sel) <- resOut
    }
    ptn@analysis@GO <- GOout
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
