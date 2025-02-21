\name{slopeFilt}
\alias{slopeFilt}
\title{Filtering \pkg{anota2seq} analysis results according to slopes for use in GSEA, GAGE, and GO term analyses
}

\description{When performing GSEA, GAGE, or GO term analysis using the output of an \pkg{anota2seq} analysis, it is often necessary to filter the input genes and log2 fold changes for the "translation" and "buffering" regulatory modes prior to performing enrichment analyses. The \code{slopeFilt} function identifies the genes with unrealistic models of changes in translation efficiency according to the slopes fitted in the per-gene \pkg{anota2seq} analysis of partial variance (APV) models. These genes are then filtered out of the analysis by providing the output of \code{slopeFilt} to the \code{genesSlopeFiltOut} argument of the \code{\link[postNet]{gseaAnalysis}}, \code{\link[postNet]{gseaPlot}}, 
\code{\link[postNet]{gageAnalysis}}, and \code{\link[postNet]{goAnalysis}} functions. See the \pkg{anota2seq} vignette for additional details on slope filtering.      
}

\usage{
slopeFilt(ads,
          regulationGen,
          contrastSel,
          minSlope,
          maxSlope)
}

\arguments{
  \item{ads}{An S4 object of class \code{Anota2seqDataSet} (resulting from the \code{\link[anota2seq]{anota2seqRun}}
  function in the \pkg{anota2seq} package). 
}
  \item{regulationGen}{The regulation effect measurement that will be taken from \code{ads}. This parameter can be
  either "translation", or "buffering". Note that slope filtering is not relevant for changes in mRNA abundance. 
}
  \item{contrastSel}{The contrast in \code{ads} where the regulation effect measurement corresponding to 
  \code{regulationGen} should be extracted from. This must be a single numeric value.
}
  \item{minSlope}{The minimum threshold to filter gene identifiers whose calculated slopes in anota2seq APV models are
  too small to be realistic. If \code{regulationGen = "translation"}, the recommended threshold is -1 (i.e., excludes
  gene identifiers with a slope < (-1)). If \code{regulationGen = "buffering"}, the recommended threshold is -2 (i.e.,
  excludes gene identifiers with a slope < (-2)). 
}
  \item{maxSlope}{The maximum threshold to filter gene identifiers whose calculated slopes in anota2seq APV models are
  too large to be realistic. If \code{regulationGen = "translation"}, the recommended threshold is 2 (i.e., excludes
  gene identifiers with a slope > 2). If \code{regulationGen = "buffering"}, the recommended threshold is 1 (i.e.,
  excludes gene identifiers with a slope > 1). 
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{A character vector of the gene identifiers that will be excluded from downstream enrichment analyses based on the selected slope filtering thresholds.
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
function (ads, regulationGen, contrastSel, minSlope, maxSlope) 
{
    check_ads(ads)
    if (!regulationGen \%in\% c("translation", "buffering")) {
        stop("For filtering slopes, 'regulationGen' should be a character vector chosen from 'translation' or 'buffering' ")
    }
    if (!is.numeric(contrastSel) | !contrastSel \%in\% seq(1, ncol(ads@contrasts), 
        1)) {
        stop("'contrastSel' should be a numeric for the desired contrast in anota2seq object")
    }
    checkSlopes(minSlope, maxSlope)
    tmpAds <- anota2seq::anota2seqGetOutput(ads, analysis = regulationGen, 
        output = "full", selContrast = contrastSel, getRVM = TRUE)
    tmpAds_slopeFilt <- tmpAds[which(tmpAds[, 1] < minSlope | 
        tmpAds[, 1] > maxSlope), ]
    genesOut <- as.character(row.names(tmpAds_slopeFilt))
    return(genesOut)
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
