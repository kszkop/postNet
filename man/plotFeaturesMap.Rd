\name{plotFeaturesMap}
\alias{plotFeaturesMap}
\title{Explore mRNA features and regulatory effects using UMAP visualizations 
}

\description{The \code{plotFeaturesMap} function implements \code{\link[umap]{umap}} to visualize the regulatory effect and features across genes.
}

\usage{
plotFeaturesMap(ptn, 
                regOnly = TRUE, 
                comparisons = NULL, 
                featSel, 
                remBinary = TRUE, 
                featCol = NULL, 
                scaled = FALSE, 
                centered = TRUE, 
                remExtreme = NULL, 
                pdfName = NULL)
}

\arguments{
  \item{ptn}{A \linkS4class{postNetData} object.
}
  \item{regOnly}{Logical specifying whether UMAPs should include only regulated genes (i.e., those included in \code{geneList}), or with all genes in the dataset. If \code{regOnly = TRUE}, only regulated genes will be used, if \code{regOnly = FALSE}, all genes will be considered. Default is \code{TRUE}. 
}
  \item{comparisons}{If \code{regOnly = TRUE}, provide a list of numeric vectors specifying pairwise comparisons between gene sets defined in the \code{regulation} or \code{geneList} slots of the \code{postNetData} object. Use 0 to denote the background set. For example, \code{list(c(0,2), c(1,2))} compares gene set 2 to the background and gene set 1 to gene set 2. The default is \code{NULL}. If both  \code{regulation} and \code{geneList} are provided, \code{geneList} is appended to the end so numbering will continue sequentially.
}
  \item{featSel}{A character vector specifying the features used in \code{\link[postNet]{featureIntegration}} to be used to generate UMAP visualizations. It is often advisable to start by exploring the data using the selected features from stepwise regression or random forest modelling.
}
  \item{remBinary}{A logical specifying if binary features should be removed when generating UMAPs. If \code{TRUE}, all binary features will be removed before generating UMAPs. If \code{FALSE}, binary features will be included. Default is \code{TRUE}.
}
  \item{featCol}{A character vector specifying which features should be coloured on the UMAPs, visualizing the relationship between the regulatory effect and the specified features. Default is \code{NULL}.
}
  \item{scaled}{A logical specifying whether the input data should be scaled prior to generating UMAPs. Default is \code{FALSE}.
}
  \item{centered}{A logical specifying whether the input data should be centered prior to generating UMAPs. Default is \code{TRUE}.
}
  \item{remExtreme}{A numeric value between \code{0} and \code{1} specifying the threshold to define extreme values to be removed from inputs prior to generating UMAPs. For example, \code{remExtreme = 0.1} removes the top 10\% of extreme values. Default is \code{NULL}.
}
  \item{pdfName}{Name to be appended to output PDF files.
}
}

\details{UMAPs are highly sensitive to input parameters. Therefore, the default settings are a suggested starting point, but may not be appropriate for all inputs. These visualizations are intended to be used as a tool for data exploration and understanding relationships between different features and regulatory effects, particularly in cases where multiple features are co-occurring in the same mRNA. 
}

\value{A PDF file with a UMAP visualization will be output for each feature specified with the \code{featCol} parameter. In each plot, the panel on the left is colored according to the regulatory effect superimposed on the UMAP generated using the features specified in the \code{featSel} parameter, while the panel on the right is coloured according to the individual feature superimposed on the same UMAP. If \code{featCol = NULL}, a PDF will be output for each feature in \code{featSel}. 

The \code{plotFeaturesMap} function also returns an updated \code{postNetData} object with a data.frame of UMAP coordinates stored in the \code{featuresMap} slot. 
}

\references{

}

\seealso{
\code{\link[umap]{umap}} \cr
\code{\link[postNet]{featureIntegration}}
}

\examples{
}
