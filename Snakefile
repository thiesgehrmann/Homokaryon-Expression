
import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
__PC_DIR__ = "%s/pipeline_components" % __INSTALL_DIR__

###############################################################################

import json
dconfig = json.load(open("%s/defaults.json"% __PC_DIR__, "r"))
dconfig.update(config)

###############################################################################

__RUN_DIR__         = os.path.abspath(dconfig["outdir"]) + "/run"
__ARA_OUTDIR__      = "%s/ara" % __RUN_DIR__
__TRANS_OUTDIR__    = "%s/trans" %  __RUN_DIR__
__MAPPING_OUTDIR__  = "%s/mapping" % __RUN_DIR__
__MARKER_OUTDIR__   = "%s/markers" % __RUN_DIR__
__QUANT_OUTDIR__    = "%s/quant" % __RUN_DIR__
__DESEQ_OUTDIR__    = "%s/deseq" % __RUN_DIR__
__KSE_OUTDIR__      = "%s/kse"% __RUN_DIR__

###############################################################################

  # Download and build the ARA toolkit
rule get_ara:
  output:
    ara = "%s/README.md" % __ARA_OUTDIR__
  conda: "%s/snakemake.yaml" % __PC_DIR__
  params:
    ara_dir = __ARA_OUTDIR__
  shell: """
    git clone https://github.com/AbeelLab/ara.git "{params.ara_dir}"
  """

rule build_ara:
  input:
    ara = rules.get_ara.output.ara
  output:
    ara = "%s/build/ara.jar" %  __ARA_OUTDIR__
  shell: """
    cd `dirname "{input.ara}"` && ant
    cd `dirname "{input.ara}"` && unzip -d build ara-development.zip
  """

###############################################################################

  # Given the FASTA and GFF files, generate the transcript sequences
rule trans:
  input:
    genome = lambda wildcards: dconfig["genomes"][wildcards.genome]["genome"],
    gff    = lambda wildcards: dconfig["genomes"][wildcards.genome]["gff"]
  output:
    fa = "%s/trans.{genome}.fa" % __TRANS_OUTDIR__
  conda: "%s/snakemake.yaml" % __PC_DIR__
  shell: """
    gffread -w "{output.fa}.orig" -g "{input.genome}" "{input.gff}"
    awk '{{ if (substr($0,1,1) == ">") {{ split($0,a," "); print a[1]}} else {{ print $0 }}}}' "{output.fa}.orig" > "{output.fa}"
  """

###############################################################################

  # Generate blast databases for each genome
rule mappingBlastDB:
  input:
    fa = lambda wildcards: "%s/trans.%s.fa" % (__TRANS_OUTDIR__, wildcards.genome)
  output:
    db = "%s/blastdb.{genome}.db"% __MAPPING_OUTDIR__
  conda: "%s/snakemake.yaml"% __PC_DIR__
  shell: """
    makeblastdb -in "{input.fa}" -dbtype nucl -out "{output.db}"
    touch "{output.db}"
  """

  # Blast genome_1 vs genome_2 and vice-versa
rule mappingBlastQuery:
  input:
    db    = lambda wildcards: expand("%s/blastdb.{genome}.db" % (__MAPPING_OUTDIR__), genome=[wildcards.genome_1]),
    trans = lambda wildcards: expand("%s/trans.{genome}.fa" % (__TRANS_OUTDIR__), genome=[wildcards.genome_2])
  output:
    res = "%s/result.{genome_1}.{genome_2}.tsv"% __MAPPING_OUTDIR__
  conda: "%s/snakemake.yaml"% __PC_DIR__
  params:
    blast_fields = dconfig["blast_fields"]
  threads: 5
  shell: """
    blastn -num_threads "{threads}" -outfmt "6 {params.blast_fields}" -evalue 0.005 -query "{input.trans}" -db "{input.db}" -out "{output.res}"
  """

  # Perform a reciprocal best blast hit to identify karyollele pairs
rule mapping:
  input:
    g1v2 = expand("%s/result.{genome_2}.{genome_1}.tsv" % __MAPPING_OUTDIR__, genome_1=["genome_1"], genome_2=["genome_2"]),
    g2v1 = expand("%s/result.{genome_1}.{genome_2}.tsv" % __MAPPING_OUTDIR__, genome_1=["genome_1"], genome_2=["genome_2"]),
    g1fa = expand("%s/trans.{genome}.fa" % __TRANS_OUTDIR__, genome=["genome_1"]),
    g2fa = expand("%s/trans.{genome}.fa" % __TRANS_OUTDIR__, genome=["genome_2"])
  output:
    mapping1 = "%s/mapping.1.fa" % __MAPPING_OUTDIR__,
    mapping2 = "%s/mapping.2.fa" % __MAPPING_OUTDIR__,
    mapping  = "%s/mapping.tsv" % __MAPPING_OUTDIR__
  params:
    blast_fields = dconfig["blast_fields"]
  run:
    import pipeline_components.utils as utils

      # Load the fasta sequences
    g1fa = utils.loadFasta(input.g1fa[0])
    g2fa = utils.loadFasta(input.g2fa[0])

      # Load the blast hits
    g12Hits = utils.indexListBy(utils.readBlastFile(input.g1v2[0], params.blast_fields), lambda x: x.qseqid)
    g21Hits = utils.indexListBy(utils.readBlastFile(input.g2v1[0], params.blast_fields), lambda x: x.qseqid)

      # For each gene, select the best hit in the other genome
    bestg12Hits = dict([ (qseqid, max(L, key=lambda x: x.bitscore)) for (qseqid,L) in g12Hits.items() ])
    bestg21Hits = dict([ (qseqid, max(L, key=lambda x: x.bitscore)) for (qseqid,L) in g21Hits.items() ])

      # Match the best hits to identify reciprocal best hits, and write the mappings to file
    matches = []
    with open(output.mapping, "w") as ofd:
      for key in bestg12Hits.keys():
        besthit = bestg12Hits[key].sseqid
        if (besthit in bestg21Hits) and (key == bestg21Hits[besthit].sseqid):
          matches.append( ((key, g1fa[key]), (besthit, g2fa[besthit])) )
          ofd.write("%s\t%s\n" % (key, besthit))
        #fi
      #efor
    #ewith

      # Write the fasta files
    utils.writeFasta([ m[0] for m in matches], output.mapping1)
    utils.writeFasta([ m[1] for m in matches], output.mapping2)


###############################################################################

  # Identify markers in the karyollele pairs
rule markers:
  input:
    ara      = rules.build_ara.output.ara,
    mapping1 = rules.mapping.output.mapping1,
    mapping2 = rules.mapping.output.mapping2,
    genomes  = [ dconfig["genomes"][genome]["genome"] for genome in dconfig["genomes"].keys() ]
  output:
    markers = "%s/markers.fa" % __MARKER_OUTDIR__
  params:
    rule_outdir = __MARKER_OUTDIR__,
    genomes = ','.join([ dconfig["genomes"][genome]["genome"] for genome in dconfig["genomes"].keys() ])
  shell: """
    java -Xmx10G -jar "{input.ara}" unique-markers "{input.mapping1},{input.mapping2}" "{params.genomes}" 21 {params.rule_outdir}/markers
    mv "{params.rule_outdir}/markers_all.fasta" {output.markers}
  """

###############################################################################

  # For each sample, quantify the markers
rule quantification:
  input:
    ara     = rules.build_ara.output.ara,
    markers = rules.markers.output.markers,
    bam     = lambda wildcards: dconfig["data"][wildcards.sample]["bam"]
  output:
    quant = "%s/quant.{sample}.tsv" % __QUANT_OUTDIR__
  params:
    paired = lambda wildcards: "--paired" if dconfig["data"][wildcards.sample].get("paired") else ""
  shell: """
    java -Xmx10G -jar "{input.ara}" snp-typer --marker "{input.markers}" -t 1 {params.paired} -o "{output.quant}" "{input.bam}"
    cat "{output.quant}.ara" \
     | grep -v '^#' \
     | cut -f1,2 \
     >  "{output.quant}"
  """

  # Summarize the markers per gene by averaging across all the markers
rule summarizeQuantPerGene:
  input:
    quant = expand("%s/quant.{sample}.tsv" % __QUANT_OUTDIR__, sample=dconfig["data"].keys())
  output:
    quant = "%s/summary.tsv" % __QUANT_OUTDIR__
  params:
    quants = [ (sample, "%s/quant.%s.tsv"% (__QUANT_OUTDIR__, sample)) for sample in dconfig["data"].keys() ]
  run:
      # Read each individual sample file
    import csv
    SQ = { sample : {} for sample in dconfig["data"].keys()  }
    for (sample, quantFile) in params.quants:
      with open(quantFile, "r") as ifd:
        reader = csv.reader(ifd, delimiter="\t")
        for row in reader:
          if len(row) != 2:
            continue
          #fi
          (marker, count) = row
          geneName = '_'.join(marker.split('_')[:-1])
          if geneName not in SQ[sample]:
            SQ[sample][geneName] = []
          #fi
          SQ[sample][geneName] = SQ[sample][geneName] + [int(count)]
        #efor
      #ewith
    #efor

      # For each gene, calculate the average counts across all the markers identified
    T = []
    samples = sorted(dconfig["data"].keys())
    for gene in SQ[samples[0]].keys():
      numMarkers = max([ len(SQ[sample][gene]) for sample in samples])
      T.append( (gene, numMarkers) + tuple([ float(sum(SQ[sample][gene]))/float(numMarkers) for sample in samples]))
    #efor

      # Write those averages to file
    with open(output.quant, "w") as ofd:
      ofd.write("#gene\tnumMarkers\t%s\n" % '\t'.join(samples))
      for geneQuant in T:
        ofd.write("%s\n" % '\t'.join([ str(x) for x in geneQuant]))
      #efor
   #ewith

###############################################################################

  # Prepare the counts file for DESeq
rule deseqInput:
  input:
    quant   = rules.summarizeQuantPerGene.output.quant,
    mapping = rules.mapping.output.mapping
  output:
    formatted = "%s/deseqInput.tsv" % __DESEQ_OUTDIR__
  params:
    minMarker = dconfig["minMarkers"]
  run:
    import csv
    M = []
    with open(input.mapping, "r") as ifd:
      reader = csv.reader(ifd, delimiter="\t")
      for (g1, g2) in reader:
        M.append((g1,g2))
      #efor
    #ewith

    GQ = {}
    with open(input.quant, "r") as ifd:
      reader = csv.reader(ifd, delimiter="\t")
      fields = []
      for row in reader:
        if row[0][0] == '#':
          fields = row[1:]
          continue
        #fi
        GQ[row[0]] = dict(zip(fields, row[1:]))
      #efor
    #ewith

    deseqRows = []
    samples = fields[1:]
    for (genome1Gene, genome2Gene) in M:
      if not((int(GQ.get(genome1Gene, { "numMarkers" : 0 })["numMarkers"]) >= params.minMarker) and (int(GQ.get(genome2Gene, { "numMarkers" : 0 })["numMarkers"]) >= params.minMarker)):
        continue
      #fi
      deseqRows.append( (genome1Gene, genome2Gene, [ (GQ[genome1Gene][sample], GQ[genome2Gene][sample]) for sample in samples ]))
    #efor

    with open(output.formatted, "w") as ofd:
      ofd.write("genegroup\t%s\n" % '\t'.join([ "%s|%s" % (genome, dconfig["data"][sample]["condition_label"])  for sample in samples for genome in [ "genome1", "genome2"]]))
      for (genome1Gene, genome2Gene, counts) in deseqRows:
        ofd.write("%s,%s\t%s\n" % (genome1Gene, genome2Gene, '\t'.join([ "%s\t%s"% (c1, c2) for (c1,c2) in counts])))
      #efor
    #ewith

  # Run DESeq
rule deseq:
  input:
    formattedQuant = rules.deseqInput.output.formatted 
  output:
    tests         = "%s/output.tests.tsv"% __DESEQ_OUTDIR__,
    counts        = "%s/output.output.tsv"% __DESEQ_OUTDIR__
  conda: "%s/snakemake.yaml" % __PC_DIR__
  params:
    deseqWrapper = "%s/deseq.R" % __PC_DIR__,
    outputPrefix = "%s/output" % __DESEQ_OUTDIR__
  shell: """
    Rscript {params.deseqWrapper} {input.formattedQuant} {params.outputPrefix}
  """

  # Using the output from DESeq, produce a table with the normalized expression for each gene in each condition in a table
rule deseqNorm:
  input:
    tests = rules.deseq.output.tests
  output:
    quant = "%s/normalizedQuant.tsv" % __DESEQ_OUTDIR__
  run:
    import csv
    Q = {}
    C = set([])
    with open(input.tests, "r") as ifd:
      reader = csv.reader(ifd, delimiter="\t")
      for row in reader:
        (condition, condG1, condG2, karyollele, basemean, basemeanG1, basemeanG2, foldchange, log2foldchange, pval, padj) = row
        karyolleles = karyollele.split(",")
        C.add(condition)
        for (k,e) in zip(karyolleles, [float(basemeanG1), float(basemeanG2)]):
          if k not in Q:
            Q[k] = {}
          #fi
          Q[k][condition] = e
        #efor
      #efor
    #ewith

    conditions = sorted(list(C))
    with open(output.quant, "w") as ofd:
      ofd.write("#gene\t%s\n" % '\t'.join(conditions))
      for gene in Q.keys():
        expr = Q.get(gene, { c:0 for c in conditions})
        ofd.write("%s\t%s\n" % (gene, '\t'.join([ str(expr.get(c, "0")) for c in conditions])))
      #efor
    #ewith

  # And provide a little function to read this file
def readDeseqNorm(quant, conditions):
  import csv
  Q = {}
  with open(quant, "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t")
    for row in reader:
      if row[0][0] == "#":
        continue
      #fi
      Q[row[0]] = dict(zip(conditions, row[1:]))
    #efor
  #ewith
  return Q
#edef

###############################################################################

  # Determine the gene Read Ratios (as described in the paper)
rule geneReadRatios:
  input:
    mapping = rules.mapping.output.mapping,
    quant   = rules.deseqNorm.output.quant
  output:
    grr = "%s/grr.tsv" % __KSE_OUTDIR__
  run:
    import pipeline_components.utils as utils
    M = utils.readMapping(input.mapping)
    C = sorted(set([config["data"][x]["condition_label"] for x in config["data"].keys()]))
    Q = readDeseqNorm(input.quant, C)

    GRR = {}
    for (g1,g2) in M:
      if (g1 not in Q) or (g2 not in Q):
        continue
      #fi
      geneGRR = []
      for c in C:
        ec1 = float(Q[g1][c]) + dconfig["pseudocount"]
        ec2 = float(Q[g2][c]) + dconfig["pseudocount"]
        geneGRR.append(ec1 / ec2)
      #efor
      GRR[(g1,g2)] = geneGRR
    #efor

    with open(output.grr, "w") as ofd:
      for (g1,g2) in GRR.keys():
        ofd.write("%s\t%s\t%s\n" % (g1,g2, '\t'.join( str(x) for x in GRR[(g1,g2)])))
      #efor
    #ewith

  # provide a function to read this file
def readGeneReadRatios(grrFile, conditions):
  import csv
  grr = {}
  with open(grrFile, "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t")
    for row in reader:
      if row[0][0] == '#':
        continue
      #fi
      grr[(row[0], row[1])] = dict(zip(conditions, [ float(x) for x in row[2:]]))
    #efor
  #ewith
  return grr
#edef

  # Determine the chromosome read ratios (as described in the paper)
rule chromosomeReadRatios:
  input:
    mapping = rules.mapping.output.mapping,
    quant   = rules.deseqNorm.output.quant,
    gff     = dconfig["genomes"]["genome_1"]["gff"]
  output:
    crr = "%s/crr.tsv" % __KSE_OUTDIR__
  run:
    import pipeline_components.utils as utils
    M = utils.readMapping(input.mapping)
    C = sorted(set([config["data"][x]["condition_label"] for x in config["data"].keys()]))
    Q = readDeseqNorm(input.quant, C)
    G = utils.readGFF3File(input.gff)
    geneChromosomes = dict([ (e.attr["ID"], e.seqid) for e in G.entries if e.type.lower() == 'mrna' ])
    

    with open(output.crr, "w") as ofd:
      ofd.write("#chromosome\t%s\n" % '\t'.join(C))
      for chromosome in G.seqids:
        crr = []
        for c in C:
          ec1 = sum([float(Q[g1][c]) for (g1,g2) in M if (g1 in Q and g2 in Q) and geneChromosomes[g1] == chromosome])
          ec2 = sum([float(Q[g2][c]) for (g1,g2) in M if (g1 in Q and g2 in Q) and geneChromosomes[g1] == chromosome])
          crr.append((ec1 + dconfig["pseudocount"]) / (ec2 + dconfig['pseudocount']))
      #efor
      ofd.write("%s\t%s\n" % (chromosome, '\t'.join([ str(x) for x in crr])))
    #efor

  # Determine the nuclear read ratios (as described in the paper)
rule nuclearReadRatios:
  input:
    mapping = rules.mapping.output.mapping,
    quant   = rules.deseqNorm.output.quant
  output:
    nrr = "%s/nrr.tsv" % __KSE_OUTDIR__
  run:
    import pipeline_components.utils as utils
    M = utils.readMapping(input.mapping)
    C = sorted(set([config["data"][x]["condition_label"] for x in config["data"].keys()]))
    Q = readDeseqNorm(input.quant, C)

    with open(output.nrr, "w") as ofd:
      ofd.write("#condition\tnrr\n")
      for c in C:
        eg1 = sum([float(Q[g1][c]) for (g1,g2) in M if (g1 in Q and g2 in Q)])
        eg2 = sum([float(Q[g2][c]) for (g1,g2) in M if (g1 in Q and g2 in Q)])
        ofd.write("%s\t%f\n" % (c, eg1/eg2))
      #efor
    #ewith

  # A wrapper rule to generate all read ratios
rule readRatios:
  input:
    grr = rules.geneReadRatios.output.grr,
    crr = rules.chromosomeReadRatios.output.crr,
    nrr = rules.nuclearReadRatios.output.nrr

###############################################################################

  # Determine the chromosome gene ratios (as described in the paper)
rule chromosomeGeneRatios:
  input:
    mapping = rules.mapping.output.mapping,
    gff     = dconfig["genomes"]["genome_1"]["gff"],
    grr     = rules.geneReadRatios.output.grr
  output:
    cgr = "%s/cgr.tsv" % __KSE_OUTDIR__
  run:
    import pipeline_components.utils as utils
    import math
    M = utils.readMapping(input.mapping)
    C = sorted(set([config["data"][x]["condition_label"] for x in config["data"].keys()]))
    GRR = readGeneReadRatios(input.grr, C)
    G = utils.readGFF3File(input.gff)
    geneChromosomes = dict([ (e.attr["ID"], e.seqid) for e in G.entries if e.type.lower() == 'mrna' ])

    with open(output.cgr, "w") as ofd:
      ofd.write("#chromosome\t%s\n" % '\t'.join(C))
      for chromosome in G.seqids:
        chrGRR = [ GRR[(g1,g2)] for (g1,g2) in M if ((g1,g2) in GRR and geneChromosomes[g1] == chromosome)]
        ngenes = len(chrGRR)
        chrCGR = []
        for c in C:
          summed = sum([math.log(q[c]) for q in chrGRR])
          chrCGR.append(math.exp( (1/max(1, float(ngenes))) * summed))
        #efor
        ofd.write("%s\t%s\n" % (chromosome, '\t'.join([ str(x) for x in chrCGR])))
     #ewith

  # Provide a function to read the output file of the previous rule
def readChromosomeGeneRatios(cgrFile, conditions):
  import csv
  cgr = {}
  with open(cgrFile, "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t")
    for row in reader:
      if row[0][0] == '#':
        continue
      #fi
      cgr[row[0]] = dict(zip(conditions, [ float(x) for x in row[1:]]))
    #efor
  #ewith
  return cgr
#edef

  # Determine the Nuclear Gene Ratios (as described in the paper)
rule nuclearGeneRatios:
  input:
    mapping = rules.mapping.output.mapping,
    gff     = dconfig["genomes"]["genome_1"]["gff"],
    cgr     = rules.chromosomeGeneRatios.output.cgr
  output:
    ngr = "%s/ngr.tsv" % __KSE_OUTDIR__
  run:
    import pipeline_components.utils as utils
    import math
    M = utils.readMapping(input.mapping)
    C = sorted(set([config["data"][x]["condition_label"] for x in config["data"].keys()]))
    CGR = readChromosomeGeneRatios(input.cgr, C)

    with open(output.ngr, "w") as ofd:
      ofd.write("#condition\tngr\n")
      for cond in C:
        nchrs = len(CGR.keys())
        summed = sum([math.log(CGR[chrom][cond]) for chrom in CGR.keys()])
        ofd.write("%s\t%f\n" % (cond, math.exp( (1/float(nchrs)) * summed)))
      #efor
    #ewith

  # A wrapper rule to generate all gene ratios.
rule geneRatios:
  input:
    cgr = rules.chromosomeGeneRatios.output.cgr,
    ngr = rules.nuclearGeneRatios.output.ngr

###############################################################################
# Generate figures... still need to do. 
