\name{BedGraphsToBigWigs}
\alias{BedGraphsToBigWigs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BedGraphsToBigWigs
}
\description{
Converts bedGraphs into bigWigs. Must write files to an output directory.
}
\usage{
BedGraphsToBigWigs(inputDir, outputDir, indexFile, bedGraphToBigWig = "bedGraphToBigWig")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputDir}{
Directory containing bedGraphs to convert to bigWigs.
}
  \item{outputDir}{
Output directory where bigWigs will be written.
}
  \item{chromSizes}{
Chromosome sizes file. Each row should be a chromosome name followed by its number of bases.
}
  \item{bedGraphToBigWig}{
Path to bedGraphToBigWig.
}
}
\details{
This function makes a system call to bedGraphToBigWig for each bedGraph file.
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
function (inputDir, outputDir, indexFile, bedGraphToBigWig = "bedGraphToBigWig")
{
    dir.create(outputDir, showWarnings = FALSE)
    for (file in list.files(path = inputDir, pattern = "*.bedGraph",
        full.names = TRUE, recursive = FALSE)) {
        outputFile <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))),
            ".bigWig")
        cmd <- paste(bedGraphToBigWig, file, indexFile, outputFile)
        system(cmd)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
