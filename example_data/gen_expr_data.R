library(polyester)
library(Biostrings)

args = commandArgs(trailingOnly=TRUE)

fasta_file = args[1]
output_dir = args[2]
stranded_arg   = args[3]
error_rate     = as.numeric(args[4])

stranded = FALSE

if (stranded_arg == "true") {
  stranded = TRUE
}

# Generate fold changes in this way
fasta = readDNAStringSet(fasta_file)
fold_change_base = c(1,1,1,1,1,1,1,1,1,1,5,4,3,2,3,1,1,1,1,1)
fold_changes = head(matrix(c(rep(fold_change_base,length(fasta)),rep(rev(fold_change_base),length(fasta))),ncol=2,),length(fasta))

readspertx = round(20 * width(fasta) / 100)

# Simulate the experiments, generating FASTA files of reads
simulate_experiment(fasta_file, strand_specific=stranded, reads_per_transcript=readspertx, 
    num_reps=c(2,2), fold_changes=fold_changes, outdir=output_dir, error_rate=error_rate)
