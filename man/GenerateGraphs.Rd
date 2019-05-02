\name{GenerateGraphs}
\alias{GenerateGraphs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GenerateGraphs
}
\description{
Generates and plots estimated partial correlation networks.
}
\usage{
GenerateGraphs(data, lambdas, rbps, approx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
An m x n matrix representing m observations (bins) and n (RBPs).
}
  \item{lambdas}{
A vector of tuning parameters to use to regularize the generated networks. length(lambdas) total networks will be generated.
}
  \item{rbps}{
A vector containing names of RBPs (in the same order as the columns are ordered). length(rbps) must equal n.
}

  \item{approx}{
Boolean, whether to use the Meinshausen and Buhlmann method to approximate the covariance matrices. Defaults to FALSE.
}
}
\details{
This function applies the nonparanormal transformation to the data before computing the correlation matrix, which it then feeds to the glassopath function. The resulting inverse covariance matrices are converted to graph weight matrices and plotted.
}
\value{
A list with the components
\item{w }{Estimated covariance matrices, has dimensions (m,n, length(lambdas))}
\item{wi }{Estimated inverse covariance matrix, has dimensions (m, n, length(lambdas))}
\item{approx}{Input argument approx}
\item{lambdas}{Vector of tuning parameters provided}
\item{errflag}{Error value, 0 indicates no error}
\item{graphs}{Weight matrices for each estimated graph, has dimensions (m, n, length(lambdas))}
}
\references{
Epskamp, S., Cramer, A.O.J., Waldorp, L.J., Schmittmann, V.D., & Borsboom, D. (2012). qgraph: Network Visualizations of Relationships in Psychometric Data. Journal of Statistical Software, 48(4). http://dx.doi.org/10.18637/jss.v048.i04.

Friedman, J., Hastie, T., & Tibshirani, R. (2007). Sparse inverse covariance estimation with the graphical lasso. Biostatistics, 0(0). 1-10.
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

set.seed(100)

x<-matrix(rnorm(50*3),ncol=3)
results <- GenerateGraphs(x, c(1e-3, 3e-3, 1e-2, 3e-2), c("RBP1", "RBP2", "RBP3"))

## The function is currently defined as
function(data, lambdas, rbps, approx=FALSE) {
	# perform non-paranormal transformation, cannot assume RBP data is normally distributed
	data.npn <- huge::huge.npn(data)

	# perform correlation
	data.cor <- qgraph::cor_auto(data.npn)

	# get graphs for each of the given lambdas
	results.path <- glasso::glassopath(data.cor, rholist=lambdas, trace=0, approx=approx)

	# set up the layout
	h <- floor(sqrt(length(lambdas)))
	w <- ceiling(length(lambdas) / h)
	layout(matrix(1:(h*w), h, w, byrow=TRUE))

	# plot each graph,
	results.path$graphs <- array(0.0, dim(results.path$wi))
	for (i in 1:length(lambdas)) {
		# convert estimated inverse-covariance matrix to weight matrix
		results.path$graphs[,,i] <- as.matrix(qgraph::wi2net(results.path$wi[,,i]))
		qgraph::qgraph(results.path$graphs[,,i], layout="spring", parallelEdge=TRUE, diag=FALSE, directed=FALSE, theme="colorblind", cut=0, title=results.path$rholist[i], labels=rbps)
	}
	# rename for clarity
	names(results.path)[4] <- "lambdas"
	results.path
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
