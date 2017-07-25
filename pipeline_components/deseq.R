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

deseq_input   = args[1];
output_prefix = args[2];

###############################################################################

countsTable <- read.delim(deseq_input, header=FALSE)
conds <- as.character(as.matrix(countsTable[1,]))[-1]
genegroups <- as.character(as.matrix(countsTable$V1))[-1]
countsTable <- countsTable[-1,-1]
countsTable <- apply(as.matrix(countsTable),2,as.numeric) # Doesn't work if there is only one gene
#countsTable <- t(as.matrix(apply(as.matrix(countsTable),2,as.numeric)))
countsTable <- round(countsTable, 0)
rownames(countsTable) <- genegroups

cond_names <- unique(sapply(conds, function(y){ strsplit(y, "|", fixed=TRUE)[[1]][[2]] }))
org_names  <- unique(sapply(conds, function(y){ strsplit(y, "|", fixed=TRUE)[[1]][[1]] }))

# https://stat.ethz.ch/pipermail/bioconductor/2010-October/035933.html
  # Remove NAs from the table
countsTable <- na.omit(countsTable)

cds <- newCountDataSet( countsTable, conds )

# Determine my own scaling factors:
# Library sizes are based on the two columns from the same sample
libSizes <- colSums(countsTable)
sampleSizes <- (libSizes[c(T,F)] + libSizes[c(F,T)])
sf <- sampleSizes / max(sampleSizes)
sizeFactors(cds) <- rep(sf, each=2)

# Estimate dispersions
# Often for smaller datasets, DESeq doesn't manage to estimate the dispersion with the normal parameters
# For example, this is the case with the example dataset.
# Therefore, we try the normal method, and if it fails, we perform a local fit, which usually works.
cds <- tryCatch({
  print("Estimating dispersions...")
  return(estimateDispersions( cds ))
}, error = function(e){
  tryCatch( {
     print("... Failed, trying with local fit...")
     return(estimateDispersions( cds ,fitType="local", locfit_extra_args=list(maxk=200)))
  }, error = function(e) {
     print("... Failed, trying with local blind fit")
     return(estimateDispersions(cds, fitType="local", method="blind"))
  })
})

# Plot the estimated dispersions
plotDispEsts(cds)
pdf(sprintf("%s.DispEsts.pdf", output_prefix))
plotDispEsts(cds)
dev.off()

Lcounts <- list()
Lres <- list()


  # For each condition
for(cond in cond_names) {
  cond_pair = paste(org_names, cond, sep='|')
    # Do the test
  res <- nbinomTest( cds, cond_pair[[1]], cond_pair[[2]] )

    # Draw an MA plot for each condition
  pdf(sprintf("%s.MAplot.%s.pdf", output_prefix, cond));
  plotMA(res)
  dev.off()

  nares <- res[(! is.na(res$padj)),]
  nA <- nrow(nares[nares$baseMeanA > nares$baseMeanB & nares$padj < 0.05, ])
  nB <- nrow(nares[nares$baseMeanA < nares$baseMeanB & nares$padj < 0.05, ])
  nE <- nrow(nares[nares$padj >= 0.05, ])

  Lcounts[[length(Lcounts)+1]] <- c(cond, paste(cond_pair, collapse='.'), nA, nB, nE)

    # Append the test output to Lres
  Lres[[length(Lres)+1]] <- cbind(condition = cond, condA = cond_pair[[1]], condB = cond_pair[[2]], res)

}

Lcounts <- t(data.frame(Lcounts))
rownames(Lcounts) <- NULL

  # Write the tests to file
write.table(Lcounts, file=sprintf("%s.output.tsv", output_prefix), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

  # Write the counts to file
write.table(do.call("rbind", Lres), file=sprintf("%s.tests.tsv", output_prefix), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE);


