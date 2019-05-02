\name{BigWigsToMatrix}
\alias{BigWigsToMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BigWigsToMatrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BigWigsToMatrix(inputDir, binsFile, bigWigAverageOverBed = "bigWigAverageOverBed", column = 4, outputDir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputDir}{
The directory containing the bigWig files to average over bins.
}
  \item{binsFile}{
The file containing the bins (observations). Each row should contain a chromosome, start, stop and a unique name. This file must be sorted by `sort -k1,1 -k2,2n`.
}
  \item{bigWigAverageOverBed}{
Path to bigWigAverageOverBed
}
  \item{column}{
Which column to use as the score for each bin. Column 3 takes the sum, 4 is the average over all bases of the bin, 5 is the average over just the covered bases. Defaults to 5 to avoid penalizing small peaks.
}
  \item{outputDir}{
Optional: Provide a directory if you would like to write the individual column files. If the directory does not exist, it will be created.
}
}
\details{
This function will make a system call to bigWigAverageOverBed for each file in the inputDir. For best performance with bigWigAverageOverBed, increase available memory.
}
\value{
BigWigsToMatrix returns a data.frame where each column corresponds to a bigWig file, and each row is a bin. The values will be sums or means depending on the bigWigAverageOverBed column selected.
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
function (inputDir, binsFile, bigWigAverageOverBed = "bigWigAverageOverBed",
    column = 4, outputDir = NULL)
{
    if (!is.null(outputDir)) {
        dir.create(outputDir, showWarnings = FALSE)
    }
    else {
        outputDir <- tempdir()
    }
    mat <- NULL
    for (file in list.files(path = inputDir, pattern = "*.bigWig",
        full.names = TRUE, recursive = FALSE)) {
        outputFile <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))),
            "_averaged.bed")
        cmd <- paste(bigWigAverageOverBed, file, binsFile, outputFile)
        system(cmd)
        result <- data.table::fread(outputFile, data.table = FALSE)
        if (!is.null(outputDir)) {
            outputFile <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))),
                ".bed")
            write.table(result[, column], outputFile, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
        if (is.null(mat)) {
            mat <- result[, column]
        }
        else {
            mat <- cbind(mat, result[, column])
        }
    }
    mat
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
