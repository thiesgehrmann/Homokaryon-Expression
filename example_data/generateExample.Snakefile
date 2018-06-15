import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

__OUTDIR__ = __INSTALL_DIR__

  # Generate the transcript sequences from the genome data files
rule genTrans:
  input:
    gff = lambda wildcards: config["genomes"][wildcards.genome]["gff"],
    genome = lambda wildcards: config["genomes"][wildcards.genome]["genome"]
  output:
    fa = "%s/trans/trans.{genome}.fa" % __OUTDIR__
  conda: "generateExample.conda.yaml"
  shell: """
    gffread -w {output.fa} -g {input.genome} {input.gff}
  """

  # combine the transcript sequences
rule transCombined:
  input:
    fa = expand("%s/trans/trans.{genome}.fa" % __OUTDIR__, genome=config["genomes"].keys())
  output:
    fa = "%s/trans/trans_combined.fa" % __OUTDIR__
  shell: """
    cat {input.fa} > {output.fa}
  """

  # Generate the fastq data using polyester, and a wrapper fasta script
rule genFastq:
  input:
    fa = rules.transCombined.output.fa
  output:
    fastq = [ ("%s/fastq/sample_01_1.fastq" % __OUTDIR__, "%s/fastq/sample_01_2.fastq" % __OUTDIR__),
              ("%s/fastq/sample_02_1.fastq" % __OUTDIR__, "%s/fastq/sample_02_2.fastq" % __OUTDIR__),
              ("%s/fastq/sample_03_1.fastq" % __OUTDIR__, "%s/fastq/sample_03_2.fastq" % __OUTDIR__),
              ("%s/fastq/sample_04_1.fastq" % __OUTDIR__, "%s/fastq/sample_04_2.fastq" % __OUTDIR__)]
  conda: "generateExample.conda.yaml"
  params:
    dir    = __INSTALL_DIR__,
    outdir = __OUTDIR__
  shell: """
    {params.dir}/fasta2fastq.sh "{input.fa}" {params.outdir}/fastq false 0.05
  """

  # Generate unaligned BAM files from the fastq files using picard and samtools
rule genSingleBam:
  input:
    fq1 = lambda wildcards: "%s/fastq/sample_%s_1.fastq" % (__OUTDIR__, wildcards.sample),
    fq2 = lambda wildcards: "%s/fastq/sample_%s_2.fastq" % (__OUTDIR__, wildcards.sample)
  output:
    bam = "%s/bam/sample_{sample}.bam" % __OUTDIR__
  conda: "generateExample.conda.yaml"
  shell: """
    picard FastqToSam F1={input.fq1} F2={input.fq2} O={output.bam}.sam SM={wildcards.sample}
    samtools sort -n -T /tmp/bam.sort.{wildcards.sample}.nnnn.bam -O bam {output.bam}.sam -o {output.bam}
  """

  # Do this for all the samples
rule genBam:
  input:
    bam = expand("%s/bam/sample_{sample}.bam" % __OUTDIR__, sample=["01", "02", "03", "04"])
