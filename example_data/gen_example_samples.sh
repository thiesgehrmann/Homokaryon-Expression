#!/usr/bin/bash

# An EXTREMELY simple example generation script.
# NOT REPRESENTATIVE OF REAL DATA
# ONLY TO DEMONSTRATE EXAMPLE

###############################################################################
  # SETTINGS
export JAVA_HOME="/usr/lib/jvm/java-1.8.0-openjdk-1.8.0.65-3.b17.el7.x86_64/"
PICARD="/home/nfs/thiesgehrmann/bulk/software/picard/dist/picard.jar"

g1_file="genome1.fasta"
g2_file="genome2.fasta"

###############################################################################

function makesample() {

  seqa=$1
  seqb=$2
  na=$3
  nb=$4
  name=$5

  for i in `seq 1 $na`; do
    makeread $seqa "${name}:a" $i
  done
  for i in `seq 1 $nb`; do
    makeread $seqb "${name}:b" $i
  done

}

function makeread() {

  seq=$1
  local name=$2
  num=$3
  len=${#seq}

  echo "@${name}:${num}"
  echo $seq
  echo "+"
  head -c $len < /dev/zero | tr '\0' 'I'
  echo -en "\n"

}

function getseqn(){
  file_name=$1
  n=$2
  cat $file_name | tr '\n' ' ' | sed -e "s/>/\n>/g" | grep -v '^$' | cut -d\  -f1 --complement | tr -d ' ' | head -n $n | tail -n1
}

###############################################################################

# Get each sequence
seqa1=`getseqn $g1_file 1`
seqa2=`getseqn $g1_file 2`
seqb1=`getseqn $g2_file 1`
seqb2=`getseqn $g2_file 2`

# Construct the sample FASTQ files
makesample $seqa1 $seqb1 100 440 sample1_1 > sample1_1.fastq
makesample $seqa2 $seqb2 150 50 sample1_1 >> sample1_1.fastq

makesample $seqa1 $seqb1 120 400 sample1_2 > sample1_2.fastq
makesample $seqa2 $seqb2 140 60 sample1_2 >> sample1_2.fastq

makesample $seqa1 $seqb1 400 100 sample2_1 > sample2_1.fastq
makesample $seqa2 $seqb2 40 120 sample2_1 >> sample2_1.fastq

makesample $seqa1 $seqb1 440 120 sample2_2 > sample2_2.fastq
makesample $seqa2 $seqb2 60 100 sample2_2 >> sample2_2.fastq

# Convert to BAM
java -Xmx10G -jar $PICARD  FastqToSam "FASTQ=./sample1_1.fastq" "OUTPUT=./sample1_1.bam" "SAMPLE_NAME=sample1_1"
java -Xmx10G -jar $PICARD  FastqToSam "FASTQ=./sample1_2.fastq" "OUTPUT=./sample1_2.bam" "SAMPLE_NAME=sample1_2"
java -Xmx10G -jar $PICARD  FastqToSam "FASTQ=./sample2_1.fastq" "OUTPUT=./sample2_1.bam" "SAMPLE_NAME=sample2_1"
java -Xmx10G -jar $PICARD  FastqToSam "FASTQ=./sample2_2.fastq" "OUTPUT=./sample2_2.bam" "SAMPLE_NAME=sample2_2"


