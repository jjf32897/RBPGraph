\name{eCLIPToBedGraph}
\alias{eCLIPToBedGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
eCLIPToBedGraph
}
\description{
Converts eCLIP narrowPeak files to bedGraphs.
}
\usage{
eCLIPToBedGraph(inputDir, outputDir, column, qualityThreshold = 1000, splitStrands = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputDir}{
Directory containing eCLIP narrowPeak files.
}
  \item{outputDir}{
Directory to which bedGraph files will be written.
}
  \item{column}{
The column number containing the score desired for analysis. Defaults to 7 (fold-enrichment).
}
  \item{qualityThreshold}{
Threshold for peak quality score to be included. Defaults to 1000 (i.e., all peaks with scores lower than 1000 will be ignored).
}
  \item{splitStrands}{
Boolean describing whether to split peaks from strands into two separate files. It is generally necessary to split strands to avoid overlapping peaks in the same file, which bedGraphToBigWig cannot handle. Splitting strands will require a call to MergeStrands on the final data matrix. Defaults to TRUE.
}
}
\details{
Outputted bedGraphs are sorted according to `sort -k1,1 -k2,2n`.
}
\value{
Returns nothing
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
function (inputDir, outputDir, column, qualityThreshold = 1000,
    splitStrands = TRUE)
{
    for (file in list.files(path = inputDir, pattern = "*.bed",
        full.names = TRUE, recursive = FALSE)) {
        d <- data.table::fread(file, data.table = FALSE)
        if (splitStrands) {
            outputFile.plus <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))),
                "_plus.bedGraph")
            outputFile.minus <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))),
                "_minus.bedGraph")
            scores.plus <- d[d[, 5] >= qualityThreshold & d[,
                6] == "+", c(1, 2, 3, column)]
            scores.minus <- d[d[, 5] >= qualityThreshold & d[,
                6] == "-", c(1, 2, 3, column)]
            write.table(setorder(scores.plus), outputFile.plus,
                sep = "\t", row.names = FALSE, col.names = FALSE,
                quote = FALSE)
            write.table(setorder(scores.minus), outputFile.minus,
                sep = "\t", row.names = FALSE, col.names = FALSE,
                quote = FALSE)
        }
        else {
            outputFile <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))),
                ".bedGraph")
            scores <- d[d[, 5] >= qualityThreshold, c(1, 2, 3,
                column)]
            write.table(setorder(scores), outputFile, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
