
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
__ANALYSIS_OUTDIR__ = "%s/analysis"% __RUN_DIR__

###############################################################################

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

# Perform a reciprocal best blast hit
rule mapping:
  input:
    g1v2 = expand("%s/result.{genome_2}.{genome_1}.tsv" % __MAPPING_OUTDIR__, genome_1=["genome_1"], genome_2=["genome_2"]),
    g2v1 = expand("%s/result.{genome_1}.{genome_2}.tsv" % __MAPPING_OUTDIR__, genome_1=["genome_1"], genome_2=["genome_2"]),
    g1fa = expand("%s/trans.{genome}.fa" % __TRANS_OUTDIR__, genome=["genome_1"]),
    g2fa = expand("%s/trans.{genome}.fa" % __TRANS_OUTDIR__, genome=["genome_2"])
  output:
    mapping1 = "%s/mapping.1.fa" % __MAPPING_OUTDIR__,
    mapping2 = "%s/mapping.2.fa" % __MAPPING_OUTDIR__,
    mapping  = "%s/mapping.tsv" % __MARKER_OUTDIR__
  params:
    blast_fields = dconfig["blast_fields"]
  run:
    import pipeline_components.utils as utils

    g1fa = utils.loadFasta(input.g1fa[0])
    g2fa = utils.loadFasta(input.g2fa[0])

    print(g1fa)
    print(g2fa)

    g12Hits = utils.indexListBy(utils.readBlastFile(input.g1v2[0], params.blast_fields), lambda x: x.qseqid)
    g21Hits = utils.indexListBy(utils.readBlastFile(input.g2v1[0], params.blast_fields), lambda x: x.qseqid)

    bestg12Hits = dict([ (qseqid, max(L, key=lambda x: x.bitscore)) for (qseqid,L) in g12Hits.items() ])
    bestg21Hits = dict([ (qseqid, max(L, key=lambda x: x.bitscore)) for (qseqid,L) in g21Hits.items() ])

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

    utils.writeFasta([ m[0] for m in matches], output.mapping1)
    utils.writeFasta([ m[1] for m in matches], output.mapping2)


###############################################################################

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

rule summarizeQuantPerGene:
  input:
    quant = expand("%s/quant.{sample}.tsv" % __QUANT_OUTDIR__, sample=dconfig["data"].keys())
  output:
    quant = "%s/summary.tsv" % __QUANT_OUTDIR__
  params:
    quants = [ (sample, "%s/quant.%s.tsv"% (__QUANT_OUTDIR__, sample)) for sample in dconfig["data"].keys() ]
  run:
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

    T = []
    samples = sorted(dconfig["data"].keys())
    for gene in SQ[samples[0]].keys():
      numMarkers = max([ len(SQ[sample][gene]) for sample in samples])
      T.append( (gene, numMarkers) + tuple([ float(sum(SQ[sample][gene]))/float(numMarkers) for sample in samples]))
    #efor

    with open(output.quant, "w") as ofd:
      ofd.write("#gene\tnumMarkers\t%s\n" % '\t'.join(samples))
      for geneQuant in T:
        ofd.write("%s\n" % '\t'.join([ str(x) for x in geneQuant]))
      #efor
   #ewith

###############################################################################
# REDO THE DESEQ PART!

###############################################################################
rule karyolleleGRR:
  input:
    quant   = rules.summarizeQuantPerGene.output.quant,
    mapping = rules.mapping.output.mapping
  output:
    grr = "%s/karyolleleFoldChanges.tsv" % __ANALYSIS_OUTDIR__
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

    GRRs = []
    samples = fields[1:]
    for (gene1, gene2) in M:
      if not((int(GQ[gene1]["numMarkers"]) >= params.minMarker) and (int(GQ[gene2]["numMarkers"]) >= params.minMarker)):
        continue
      #fi
      GRRs.append( (gene1, gene2) + tuple( [ float(GQ[gene1][sample]) / float(GQ[gene2][sample]) for sample in samples ]))
    #efor

    with open(output.grr, "w") as ofd:
      ofd.write("#genome1_gene\t#genome2_gene\t%s\n" % '\t'.join(samples))
      for grr in GRRs:
        ofd.write("%s\t%s\t%s\n" % (grr[0], grr[1], '\t'.join([ str(f) for f in grr[2:]])))
      #efor
    #ewith
 
        
