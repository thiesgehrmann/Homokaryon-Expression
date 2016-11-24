#!/bin/sh

###############################################################################
# CONFIGURATION

ARALOC="/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/show_marcel/ara/build/ara.jar"

###############################################################################

function usage(){

  echo "$0 Usage"
  echo "  $0 <genome_files> <gff_files> <output_dir> <data_file> <paired_end>"
  echo "";
  echo "genome_files:  Comma-separated list of FASTA files, representing each homokaryon/allele";
  echo "gff_files:     Comma-separated list of GFF files, representing each homokaryon/allele";
  echo "output_dir:    Directory to which output should be written";
  echo "data_file:     File with BAM files and samples described"
  echo "paired_end:    [Y|N] Is the paired end?"

}

###############################################################################

if [  $# -lt 4 ]; then

  usage $0;
  exit 1;
fi

genome_files="$1";
shift;
gff_files="$1";
shift;
output_dir="$1";
shift;
data_file="$1";
shift
paired="$1";
shift
tissue_figure="$1"
shift


mkdir -p $output_dir

n_genomes=`echo $genome_files | tr ',' '\n' | wc -l`

###############################################################################
# Mapping - Identify pairs
#./mapping.py $genome_files $gff_files "$output_dir"

mapping_files_produced=`for x in $(seq 0 $((n_genomes-1))); do echo "$output_dir/mapping_${x}.fasta"; done | tr '\n' ',' | sed -e 's/,$//'`

###############################################################################
# Identify unique sequences that uniquely identify genes between groups of genes
#~/scala-2.11.7/bin/scala -J-Xmx30G uniqueProbes "$mapping_files_produced" "$genome_files" 21 "$output_dir/UniqueMarkers"
#java -jar -Xmx10G ../ara/build/ara.jar unique-markers "$mapping_files_produced" "$genome_files" 21 "$output_dir/UniqueMarkers"

# Merge all these unique tags into one file
find "$output_dir"  | grep -e "UniqueMarkers_[0-9]\+[.]fasta" | xargs cat > $output_dir/UniqueMarkers_all.fasta

###############################################################################
# Determine abundance of each tag in each sample

#java -jar -Xmx10G ../macaw/ara.jar snp-typer --marker uniqueProbes_all.fasta -t 1 -o ./probeCounts.tsv /home/nfs/thiesgehrmann/bulk/unaligned_bam/agabi/YFB_F2737_R1.bam
#sauto short -cmd "java -jar -Xmx30G ../ara/build/ara.jar unique-markers /home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/mapping_0.fasta,/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/mapping_1.fasta $genome_files 21 /home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/uniqueProbes"
#for fastq_file in "$@";
#  print
#  # Run ARA to count for each tag

cat $data_file | grep -v '^$' | grep -v '^[ \t]*#' | sed -e 's/[ \t]\+/ /g' | while read data_line; do
  sample_name=`echo "$data_line" |  cut -d\  -f1`;
  replicate_id=`echo "$data_line" | cut -d\  -f2`;
  replicate_label=`echo "$data_line" | cut -d\  -f3`;
  bam_file=`echo "$data_line" | cut -d\   -f4`;
  if [ "$paired" = 'Y' ]; then
    paired='--paired'
  else
    paired='';
  fi

  echo $data_line
  echo $sample_name
  echo $replicate_id
  echo $replicate_label
  echo $bam_file
  
  java -jar -Xmx10G /home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/show_marcel/ara/build/ara.jar snp-typer --marker "$output_dir/UniqueMarkers_all.fasta" -t 1 $paired -o "$output_dir/MarkerCounts_${sample_name}.tsv" "$bam_file"
done


###############################################################################
# Analysis of tag counts

./analysis.py $data_file $genome_files $gff_files $output_dir $tissue_figure
###############################################################################

