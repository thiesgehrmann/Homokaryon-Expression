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

if [  $# -lt 3 ]; then

  usage $0;
  exit 1;
fi

genome_files="$1";
shift;
gff_files="$1";
shift;

genome_files="/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h39_1/AgabiH39_1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h97_1/agabiH97_1.assembly.fasta"
gff_files="/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h39_1/AgabiH39_1.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h97_1/AgabiH97_1.gff3"
output_dir='/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/U1'
fastq_files="$@";

genome_files="/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.assembly.fasta"
gff_files="/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.genes.gff3"
output_dir='/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/A15'
fastq_files="$@";

mkdir -p $output_dir

n_genomes=`echo $genome_files | tr ',' '\n' | wc -l`

###############################################################################
# Mapping - Identify pairs
./mapping.py $genome_files $gff_files "$output_dir"

mapping_files_produced=`for x in $(seq 0 $((n_genomes-1))); do echo "$output_dir/mapping_${x}.fasta"; done | tr '\n' ',' | sed -e 's/,$//'`

###############################################################################
# Identify unique sequences that uniquely identify genes between groups of genes
#~/scala-2.11.7/bin/scala -J-Xmx30G uniqueProbes "$mapping_files_produced" "$genome_files" 21 "$output_dir/UniqueMarkers"
java -jar -Xmx10G ../ara/build/ara.jar unique-markers "$mapping_files_produced" "$genome_files" 21 "$output_dir/UniqueMarkers"

# Merge all these unique tags into one file
find "$output_dir"  | grep -e "UniqueMarkers_[0-9]\+[.]fasta" | xargs cat > $output_dir/UniqueMarkers_all.fasta

###############################################################################
# Determine abundance of each tag in each sample

#java -jar -Xmx10G ../macaw/ara.jar snp-typer --marker uniqueProbes_all.fasta -t 1 -o ./probeCounts.tsv /home/nfs/thiesgehrmann/bulk/unaligned_bam/agabi/YFB_F2737_R1.bam
#sauto short -cmd "java -jar -Xmx30G ../ara/build/ara.jar unique-markers /home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/mapping_0.fasta,/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/mapping_1.fasta $genome_files 21 /home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/uniqueProbes"
#for fastq_file in "$@";
#  print
#  # Run ARA to count for each tag

bamdir="/home/nfs/thiesgehrmann/bulk/unaligned_bam/agabi";

#for bamfile in `ls $bamdir | grep ".bam$"`; do
#  bamname=`echo $bamfile | sed -e 's/[.]bam$//'`;
#  echo $bamdir/$bamfile, $bamname
#  samtools sort -n -o $bamdir/${bamname}.sorted.bam -@4 $bamdir/$bamfile;
#done

for bamfile in `ls $bamdir | grep ".bam$"`; do
  bamname=`echo $bamfile | sed -e 's/[.]bam$//'`;
  java -jar -Xmx10G ../ara/build/ara.jar snp-typer --marker "$output_dir/UniqueMarkers_all.fasta" -t 1 -o --paired "$outpit_dir/MarkerCounts_${bamname}.tsv" "$bamdir/$bamfile"
done


###############################################################################
# Analysis of tag counts


###############################################################################
