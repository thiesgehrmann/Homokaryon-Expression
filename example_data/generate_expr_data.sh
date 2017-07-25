#!/usr/bin/env bash

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

genes="$1"
outdir="$2"
stranded="$3"
errorrate="$4"

mkdir -p "$outdir"
Rscript $SCRIPTDIR/gen_expr_data.R $genes $outdir $stranded $errorrate

find $outdir \
 | grep -e '[.]fasta$' \
 | while read fn; do
  newname=`echo $fn | sed -e 's/[.]fasta$/.fastq/'`
  cat $fn | $SCRIPTDIR/fasta2fastq.sh > $newname
  rm $fn
done

