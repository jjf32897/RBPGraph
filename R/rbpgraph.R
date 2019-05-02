# Generates and plots estimated partial correlation networks for each provided lambda value
GenerateGraphs <- function(data, lambdas, rbps, approx=FALSE) {
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
		qgraph::qgraph(results.path$graphs[,,i], layout="spring", parallelEdge=TRUE, diag=FALSE, directed=FALSE, theme="colorblind",
			cut=0, title=results.path$rholist[i], labels=rbps)
	}
	# rename for clarity
	names(results.path)[4] <- "lambdas"
	results.path
}

# Executes bigWigAverageOverBed
BigWigsToMatrix <- function(inputDir, binsFile, bigWigAverageOverBed="bigWigAverageOverBed", column=5, outputDir=NULL) {
	# create output directory if it doesn't exist, get temporary output directory if none given
	if (!is.null(outputDir)) {
		dir.create(outputDir, showWarnings=FALSE)
	} else {
		outputDir <- tempdir()
	}

	# final matrix
	mat <- NULL

	# iterate over all bigWigs in inputDir
	for (file in list.files(path=inputDir, pattern="*.bigWig", full.names=TRUE, recursive=FALSE)) {
		# run bigWigAverageOverBed
		outputFile <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))), '_averaged.bed')
		cmd <- paste(bigWigAverageOverBed, file, binsFile, outputFile)
		system(cmd)

		# read in the matrix just created
		result <- data.table::fread(outputFile, data.table=FALSE)

		# if there was a given outputDir, write the column to a file
		if (!is.null(outputDir)) {
			outputFile <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))), '.bed')
			write.table(result[,column], outputFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
		}

		# if this is the first file, set mat to it, otherwise append new column
		if (is.null(mat)) {
			mat <- result[,column]
		} else {
			mat <- cbind(mat, result[,column])
		}
	}
	mat
}

# Converts bedGraphs into bigWigs. Writes to an output directory, returns nothing.
BedGraphsToBigWigs <- function(inputDir, outputDir, chromSizes, bedGraphToBigWig="bedGraphToBigWig") {
	# attempt to create the output directory
	dir.create(outputDir, showWarnings=FALSE)

	# iterate over all bedGraphs in inputDir
	for (file in list.files(path=inputDir, pattern="*.bedGraph", full.names=TRUE, recursive=FALSE)) {
		# run bedGraphToBigWig
		outputFile <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))), '.bigWig')
		cmd <- paste(bedGraphToBigWig, file, chromSizes, outputFile)
		system(cmd)
	}
}

# Takes eCLIP peak files (BED format) and generates bedGraph files with the peaks and their specified scores.
# Users can specify the quality threshold for the peaks, which column to use as the score, and also whether
# to split strands. Writes to an output file, returns nothing.
# Without splitting strands, you may have overlapping peaks, which is BAD!
eCLIPToBedGraph <- function(inputDir, outputDir, column=7, qualityThreshold=1000, splitStrands=TRUE) {
	for (file in list.files(path=inputDir, pattern="*.bed", full.names=TRUE, recursive=FALSE)) {
		d <- data.table::fread(file, data.table=FALSE)
		# get only the ch, start, stop, and column score for all peaks exceeding the quality threshold
		# write them to the output directory
		if (splitStrands) {
			# separate plus and minus strands
			outputFile.plus <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))), '_plus.bedGraph')
			outputFile.minus <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))), '_minus.bedGraph')
			scores.plus <- d[d[,5] >= qualityThreshold & d[,6] == '+',c(1, 2, 3, column)]
			scores.minus <- d[d[,5] >= qualityThreshold & d[,6] == '-',c(1, 2, 3, column)]
			# sort them
			write.table(setorder(scores.plus), outputFile.plus, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
			write.table(setorder(scores.minus), outputFile.minus, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
		} else {
			outputFile <- paste0(file.path(outputDir, tools::file_path_sans_ext(basename(file))), '.bedGraph')
			scores <- d[d[,5] >= qualityThreshold,c(1, 2, 3, column)]
			write.table(setorder(scores), outputFile, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
		}
	}
}

# If matrix formed from files where + and - strands were split, they need to be merged. Assuming adjacent columns
# are from the same RBP, we merge adjacent columns, either by taking the average or summing.
MergeStrands <- function(mat, fn="mean") {
	divisor <- 2
	if (fn == "sum") {
		divisor <- 1
	}
	mat <- (mat[,seq(1, ncol(mat), by=2)] + mat[, seq(2, ncol(mat), by=2)]) / divisor
	mat
}