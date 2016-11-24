#!/bin/sh

###############################################################################
#U1

genome_files="/home/nfs/thiesgehrmann/groups/w/phd/data/U1_homokaryons/AgabiH39_2.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/U1_homokaryons/AgabiH97_2.assembly.fasta"
gff_files="/home/nfs/thiesgehrmann/groups/w/phd/data/U1_homokaryons/AgabiH39_2.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/U1_homokaryons/AgabiH97_2.genes.gff3"
output_dir='/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/U1'
data_file="/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/agabi_tissue_data.tsv"
paired="Y"
tissue_figure="/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/homokaryon_activity.svg"

./pipeline.sh $genome_files $gff_files $output_dir $data_file $paired $tissue_figure

###############################################################################
# U1 COMPOST DATA

genome_files="/home/nfs/thiesgehrmann/groups/w/phd/data/U1_homokaryons/AgabiH39_2.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/U1_homokaryons/AgabiH97_2.assembly.fasta"
gff_files="/home/nfs/thiesgehrmann/groups/w/phd/data/U1_homokaryons/AgabiH39_2.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/U1_homokaryons/AgabiH97_2.genes.gff3"
output_dir='/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/U1_compost'
data_file="/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/agabi_compost_data.tsv"
paired="N"
tissue_figure="";

./pipeline.sh $genome_files $gff_files $output_dir $data_file $paired $tissue_figure


###############################################################################
# A15

genome_files="/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.assembly.fasta"
gff_files="/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.genes.gff3"
output_dir='/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/A15'
data_file="/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/agabi_tissue_data.tsv"
paired="Y"
tissue_figure="/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/homokaryon_activity.svg"

./pipeline.sh $genome_files $gff_files $output_dir $data_file $paired $tissue_figure


###############################################################################
# A15 COMPOST DATA

genome_files="/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.assembly.fasta"
gff_files="/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.genes.gff3"
output_dir='/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/A15_compost'
data_file="/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/agabi_compost_data.tsv"
paired="N"
tissue_figure="";

./pipeline.sh $genome_files $gff_files $output_dir $data_file $paired $tissue_figure


###############################################################################
# A15 COMPOST DATA (NODUP)

genome_files="/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.assembly.fasta"
gff_files="/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.genes.gff3"
output_dir='/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/A15_compost_nodup'
data_file="/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/agabi_compost_nodup.tsv"
paired="N"
tissue_figure="";

./pipeline.sh $genome_files $gff_files $output_dir $data_file $paired $tissue_figure

