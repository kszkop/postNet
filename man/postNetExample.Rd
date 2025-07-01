\name{postNetExample}
\alias{postNetExample}
\docType{data}
\title{
A Capitalized Title for the Data Set
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("postNetExample")}
\format{
  The format is:
List of 4
 $ geneList  :List of 2
  ..$ translationUp  : chr [1:1445] "DDX6" "CLK1" "DUSP1" "TXNIP" ...
  ..$ translationDown: chr [1:1378] "PPDPF" "NCKAP1" "AP2M1" "ABHD12" ...
 $ background: chr [1:9555] "A4GALT" "AAAS" "AACS" "AADAT" ...
 $ effect    : Named num [1:9555] -2.095 0.275 -0.57 -0.608 0.576 ...
  ..- attr(*, "names")= chr [1:9555] "A4GALT" "AAAS" "AACS" "AADAT" ...
 $ ads_data  :List of 3
  ..$ phenoVec: chr [1:8] "Control" "Control" "Control" "Control" ...
  ..$ dataP   : num [1:100, 1:8] 6.15 7.85 8.62 6.15 5.8 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:100] "ACOX1" "ACTR2" "AGRN" "APOBEC3C" ...
  .. .. ..$ : chr [1:8] "riboSeq_Control_rep1" "riboSeq_Control_rep2" "riboSeq_Control_rep3" "riboSeq_Control_rep4" ...
  ..$ dataT   : num [1:100, 1:8] 5.77 8.23 9.39 5.34 6.4 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:100] "ACOX1" "ACTR2" "AGRN" "APOBEC3C" ...
  .. .. ..$ : chr [1:8] "rnaSeq_Control_rep1" "rnaSeq_Control_rep2" "rnaSeq_Control_rep3" "rnaSeq_Control_rep4" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(postNetExample)
## maybe str(postNetExample) ; plot(postNetExample) ...
}
\keyword{datasets}
