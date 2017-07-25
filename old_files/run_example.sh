genome_files="./example_data/genome1.fasta,./example_data/genome2.fasta"
gff_files="./example_data/genome1.gff3,./example_data/genome2.gff3"
output_dir='./example_output'
data_file="./example_data/dataset.tsv"
paired="N"

./pipeline.sh $genome_files $gff_files $output_dir $data_file $paired
