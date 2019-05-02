\name{MergeStrands}
\alias{MergeStrands}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MergeStrands
}
\description{
Merges adjacent columns of a matrix (assumed to be opposite strands of the same RBP) by combining the values either by sum or average.
}
\usage{
MergeStrands(mat, fn = "mean")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
Final data matrix, where each column corresponds to a strand of an RBP. If this function is being called, it is implied that eCLIPToBedGraph split strands; i.e., each eCLIP file yielded two bedGraphs.
}
  \item{fn}{
Describes how to combine the scores. Default is "mean", but can be changed to "sum".
}
}
\details{
This function should be called on the final matrix if and only if splitStrands was set to TRUE in eCLIPToBedGraph. The input matrix should have 2n columns, where n is the total number of RBPs. As long as the file names created by upstream RBPGraph functions were left unchanged, opposing strands from the same RBP should be adjacent to one another, which is the assumption this function makes. Confirm that it is the case where files in the inputDir to bigWigsToMatrix were ordered such that each pair of adjacent files corresponds to the same RBP.
}
\value{
Returns an m x n matrix, where the input was m x 2n. Adjacent pairs of columns are combined either by averaging the values or summing them.
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
function (mat, fn = "mean")
{
    divisor <- 2
    if (fn == "sum") {
        divisor <- 1
    }
    mat <- (mat[, seq(1, ncol(mat), by = 2)] + mat[, seq(2, ncol(mat),
        by = 2)])/divisor
    mat
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
