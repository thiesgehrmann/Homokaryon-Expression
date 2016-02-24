# Using http://dwheelerau.com/2013/04/15/how-to-use-deseq-to-analyse-rnaseq-data/

library("DESeq")

###############################################################################

printf <- function(...) invisible(print(sprintf(...)))

plotDispEsts <- function( cds ){

  plot(rowMeans( counts( cds, normalized=TRUE ) ), fitInfo(cds)
    $perGeneDispEsts, pch = '.', log="xy", ylab="Dispersion",
    xlab="Mean of normalized counts")
  xg = 10^seq( -.5, 5, length.out=300 )
  lines( xg, fitInfo(cds)$dispFun( xg ), col="red")

}

###############################################################################

args <- commandArgs(trailingOnly = TRUE)

deseq_input = args[1];
output_dir  = args[2];

###############################################################################

countsTable <- read.delim(deseq_input, header=FALSE)
conds <- as.character(as.matrix(countsTable[1,]))[-1]
genegroups <- as.character(as.matrix(countsTable$V1))[-1]
countsTable <- countsTable[-1,-1]
countsTable <- apply(as.matrix(countsTable),2,as.numeric)
countsTable <- round(countsTable, 0)
rownames(countsTable) <- genegroups

cond_names <- unique(sapply(conds, function(y){ strsplit(y, "|", fixed=TRUE)[[1]][[2]] }))
org_names  <- unique(sapply(conds, function(y){ strsplit(y, "|", fixed=TRUE)[[1]][[1]] }))

# I may need to estimate my own size factors for these myself
# https://stat.ethz.ch/pipermail/bioconductor/2010-October/035933.html
  # Remove NAs from the table
countsTable <- na.omit(countsTable)

cds <- newCountDataSet( countsTable, conds )

# Here we observe the differences in expression between the two samples.
# We sum the total read counts for each homokaryon in each sample, and then we compare the two within each sample.
# We look at the distribution of these values for the two replicates.
# This is put in condLibRatios
libSizes <- colSums(countsTable)
sampleLibRatios <- libSizes[c(TRUE, FALSE)] / libSizes[c(FALSE, TRUE)]
condLibRatios <- data.frame(r1=sampleLibRatios[c(rep(F, length(conds)/4), rep(T, length(conds)/4))], r2=sampleLibRatios[rev(c(rep(F, length(conds)/4), rep(T, length(conds)/4)))])
condLibRatios <- transform(condLibRatios, AVG=apply(condLibRatios,1, mean, na.rm = TRUE), SD=apply(condLibRatios,1, sd, na.rm = TRUE), VAR=apply(condLibRatios,1, var, na.rm = TRUE))
rownames(condLibRatios) <- cond_names

write.table(condLibRatios, file=sprintf("%s/deseq_condlibratios.tsv", output_dir), sep='\t', row.names=TRUE, col.names=FALSE, quote=FALSE)

# Determine my own scaling factors:
# Library sizes are based on the two columns from the same sample
# And then estimate the dispersions
sampleSizes <- (libSizes[c(T,F)] + libSizes[c(F,T)])
sf <- sampleSizes / max(sampleSizes)
sizeFactors(cds) <- rep(sf, each=2)
cds <- estimateDispersions( cds )


plotDispEsts(cds)
jpeg(sprintf("%s/DispEsts_DESeq.jpg", output_dir))
plotDispEsts(cds)
dev.off()

Lcounts <- list()
Lres <- list()
for(cond in cond_names) {
  cond_pair = paste(org_names, cond, sep='|')
  res <- nbinomTest( cds, cond_pair[[1]], cond_pair[[2]] )
  pdf(sprintf("%s/MAplot_%s.pdf", output_dir, cond));
  plotMA(res)
  dev.off()
  nares <- res[(! is.na(res$padj)),]
  nA <- nrow(nares[nares$baseMeanA > nares$baseMeanB & nares$padj < 0.05, ])
  nB <- nrow(nares[nares$baseMeanA < nares$baseMeanB & nares$padj < 0.05, ])
  nE <- nrow(nares[nares$padj >= 0.05, ])

  Lcounts[[length(Lcounts)+1]] <- c(cond, paste(cond_pair, collapse='.'), nA, nB, nE)

  Lres[[length(Lres)+1]] <- cbind(condition = cond, condA = cond_pair[[1]], condB = cond_pair[[2]], res)

}
Lcounts <- t(data.frame(Lcounts))
rownames(Lcounts) <- NULL
write.table(Lcounts, file=sprintf("%s/deseq_output.tsv", output_dir), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(do.call("rbind", Lres), file=sprintf("%s/deseq_tests.tsv", output_dir), sep='\t', col.names=FALSE, row.names=TRUE, quote=FALSE);


