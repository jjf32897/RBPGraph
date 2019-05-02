# RBPGraph
An R package to estimate RBP co-binding networks from eCLIP narrowPeak data.
Provides a series of functions to automate the process from eCLIP files -> bedGraphs -> bigWigs -> RBP-bin observation matrix -> graphs.

Depends on the glasso, qgraph, and huge packages.

Please refer to the PDF manual in this repository for more details.

Here is an example of a typical workflow:

```
# define data and software paths
eCLIP.dir <- "ENCODE/eCLIP_peaks/"
bedGraph.dir <- "GGM_data/bedGraphs/"
bigWig.dir <- "GGM_data/bigWigs/"
chromSizes <- "annotations/hg38.chrom.sizes"
my.bins <- "annotations/hg38.bins.100.bed"
bg2bw <- "software/bedGraphToBigWig"
bwaob <- "software/bigWigAverageOverBed"

# define parameters for graph generation/plotting
lambdas <- c(1e-3, 3e-3, 1e-2, 3e-2)
my.rbps <- c("AGGF1", "BCCIP", "BUD13")

# call RBPGraph functions step-by-step
eCLIPToBedGraph(eCLIP.dir, bedGraph.dir)
BedGraphsToBigWigs(bedGraph.dir, bigWig.dir, chromSizes, bedGraphToBigWig=bg2bw)
d <- BigWigsToMatrix(bigWig.dir, my.bins, bigWigAverageOverBed=bwaob)
results <- GenerateGraphs(d, lambdas, my.rbps)
```

Thank you to Jing Zhang, Mark Gerstein, Jason Liu, and Donghoon Lee for their help on this project!
