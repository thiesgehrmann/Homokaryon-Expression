#!/bin/sh

function usage(){

  echo "$0 Usage"
  echo "  $0 <genome_files> <gff_files> <fastq_files_1> [fastq_files_2] ... [fastq_files_n]"
  echo "";
  echo "genome_files:  Comma-separated list of FASTA files, representing each homokaryon/allele";
  echo "gff_files:     Comma-separated list of GFF files, representing each homokaryon/allele";
  echo "fastq_files_n: FASTQ file of each sample, if paired end, then comma-separated";

}

###############################################################################

if [  $# -lt 3 ];

  usage $0;
  exit 1;
fi

genome_files="$1";
shift;
gff_files="$1";
shift;

genome_files="/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h39_1/AgabiH39_1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h97_1/agabiH97_1.assembly.fasta"
gff_files="/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h39_1/AgabiH39_1.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h97_1/AgabiH97_1.gff3"


###############################################################################
# Mapping - Identify pairs


###############################################################################
# Identify unique sequences that uniquely identify genes between groups of genes


###############################################################################
# Determine abundance of each tag in each sample
for fastq_file in "$@";
  print
  # Run ARA to count for each tag


###############################################################################
# Analysis of tag counts


###############################################################################
