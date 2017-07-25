#!/bin/env python

import sys

import numpy as np
from ibidas import *
from scipy.stats import t
import python_utils as utils;
import python_svg_utils as svg_utils;
import pandas as pd
import matplotlib
import matplotlib.pylab as plt;
import seaborn as sns;
import genome_drawing as gd;
import genome_drawing_split as gd2
import math

###############################################################################

# The structure of the data

#sys.argv = ['',
#            'agabi_tissue_data.tsv',
#            "/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.assembly.fasta",
#            "/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.genes.gff3",
#            'A15',
#            'homokaryon_activity.svg']

#sys.argv = [ '',
#             'agabi_compost_nodup.tsv',
#             "/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.assembly.fasta",
#             "/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.genes.gff3",
#             'A15_compost_nodup']



data_file  = Read(sys.argv[1]).Cast(str, int, str, str) / ('sample_name', 'replicate_id', 'replicate_label', 'bam_file');
genomes    = [ Read(x, format='fasta') for x in sys.argv[2].split(',')]
gffs       = [ Read(x, format='gff') for x in sys.argv[3].split(',')]
output_dir = sys.argv[4];
tissue_figure = sys.argv[5] if len(sys.argv) == 6 else None;

sample_names = data_file.Sort(_.replicate_id).Get(_.sample_name)()
sample_labels = data_file.Sort(_.replicate_id).Get(_.replicate_id)()
label_names   = data_file.Get(_.replicate_id, _.replicate_label).Unique().Sort(_.replicate_id).replicate_label()

###############################################################################

# Gather all the read counts, calculate average per gene, and put into a big table

allD = [];
allDs = []
for i, label_name in enumerate(label_names):
  samples    = sample_names[sample_labels == i]
  sample_ids = np.where(sample_labels == i)[0]
  print samples;
  Ds = [ Read("%s/MarkerCounts_%s.tsv.ara" % (output_dir, x), verbose=False).Detect() for x in samples];
  Ds = [ x.To(_.marker, Do=_.Each(lambda x: '_'.join(x.split('_')[:-1])).Cast(str)) for x in Ds ]
  Ds = [ x.ReplaceMissing().GroupBy(_.marker) for x in Ds]
  Ds = [ x.To(_.read_depth, Do=_.Mean()) for x in Ds]
  Ds = [ x.To(_.p_value, Do=_.Count()) for x in Ds ]
  Ds = [ x.Without(_.a_p) for x in Ds ]
  print label_name
  Ds = [ x.Get(_.marker, _.marker.Each(lambda x: label_name).Cast(str), (_.read_depth * 0  + s_id).Cast(int), _.read_depth, _.p_value).Copy() for (s_id, x) in zip(sample_ids.tolist(), Ds) ]
  Ds = [ x / ('geneid', 'sample_name', 'replicate_id', 'avg_read_depth', 'nmarkers') for x in Ds ]
  Dstacked = Ds[0]
  print Ds[0].Names
  for d in Ds[1:]:
    Dstacked = Dstacked | Stack | d
  #efor
  allD.append(Dstacked);
  allDs = allDs + Ds;
#efor

D = allD[0];
for d in allD[1:]:
  D = D | Stack | d;
#efor
D = D.Copy()

D = D.To(_.avg_read_depth, Do=_.Each(lambda x: math.ceil(x)).Cast(int)).Copy();

###############################################################################
# Examine the correlation between average expression per gene and the average cg content per gene markers

mapping   = Read('%s/mapping.tsv' % output_dir) / ('id0', 'id1')
markers_0 = (Read('%s/UniqueMarkers_0.tsv' % output_dir) / ('id0', 'loc', 'seq0')) | Match(_.id0, merge_same='equi') | mapping
markers_1 = (Read('%s/UniqueMarkers_1.tsv' % output_dir) / ('id1', 'loc', 'seq1')) | Match(_.id1, merge_same='equi') | mapping

def cg_content(str):
  l  = len(str)
  cg = len([c for c in str if c in [ 'C', 'G', 'c', 'g']])
  return float(cg) / l
#edef

G_CG_0 = (markers_0.To(_.seq0, Do=_.Each(lambda x: cg_content(x)).Cast(float))).GroupBy(_.id0)[_.seq0.Count() > 5].Get(_.id0, _.seq0.Mean(), _.id1[0])
G_CG_1 = (markers_1.To(_.seq1, Do=_.Each(lambda x: cg_content(x)).Cast(float))).GroupBy(_.id1)[_.seq1.Count() > 5].Get(_.id1, _.seq1.Mean(), _.id0[0])

M_CG = (G_CG_0 | Match(_.id0, merge_same='all', jointype='inner') | G_CG_1).Get(_.id0, _.id1, _.seq0, _.seq1)

M_CG | Match(_.id0, _.geneid, merge_same='equi', jointype='left') | (D.Get(_.geneid, _.sample_name, _.replicate_id, _.avg_read_depth))

M_CG_RD = ((M_CG | Match(_.id0, _.geneid, merge_same='all', jointype='left') | (D.Get(_.geneid, _.sample_name, _.replicate_id, _.avg_read_depth / 'rd0'))) | Match((_.id1, _.replicate_id), (_.geneid, _.replicate_id), merge_same='all', jointype='left') | (D.Get(_.geneid, _.sample_name, _.replicate_id, _.avg_read_depth / 'rd1'))).Get(_.id0, _.id1, _.seq0, _.seq1, _.rd0, _.rd1, _.sample_name, _.replicate_id).Copy()

plt.cla(); plt.clf();
C_mcgrd = M_CG_RD[_.replicate_id == 0].Get((_.seq0 - _.seq1).Cast(float), ((_.rd0 - _.rd1).Cast(float) / (_.rd0 + _.rd1)).Cast(float))()
sns.regplot(C_mcgrd[0], C_mcgrd[1], fit_reg=False)
p = np.corrcoef(C_mcgrd[0], C_mcgrd[1])[1][0]
plt.xlabel("Difference between CG content");
plt.ylabel("Fraction difference between read counts");
plt.title("Vegetative (p=%f)" % p)
#plt.axes().set_yscale("log")
#plt.axes().set_xscale("log")
plt.savefig('%s/CG_read_depth.png' % output_dir);
plt.savefig('%s/CG_read_depth.svg' % output_dir);


###############################################################################
# Based on the mapping, create a view where we see both

M = Read('%s/mapping.tsv' % (output_dir))

Mgroup = M.Get(_.Get(tuple(M.Names)).To(0, Do=_.Each(lambda x: ','.join(x)).Cast(str) / 'group'), *M.Names);
Mflat = Mgroup.Get(_.group, M.Names[0])
for orgname in M.Names[1:]:
  Mflat = Mflat | Stack | Mgroup.Get(_.group, orgname);
#efor
Mflat = Mflat / ('genegroup', 'geneid');

Dgenegroup = Mflat | Match(_.geneid, _.geneid) | D
Dgenegroup = Dgenegroup.GroupBy(_.genegroup)[_.nmarkers.Min() > 5].Flat();
Dgenegroup = Dgenegroup.Sort(_.geneid, _.replicate_id);

###############################################################################
###############################################################################
###############################################################################
  # Prepare DESEQ input

###############################################################################
# Gene level

fieldnames  = Dgenegroup.GroupBy(_.genegroup).Sort(_.geneid, _.replicate_id)[0].Get(_.geneid, _.sample_name).To(_.geneid, Do=_.Each(lambda x: x.split('|')[0]).Cast(str)).Get(_.geneid + '|' + _.sample_name / 'fieldnames')().tolist()
fieldnames  = tuple([ x.lower() for x in ['genegroup'] + fieldnames ])
DESEQ_genes = Rep([ (x,) + tuple(y.tolist()) for (x,y) in zip(*Dgenegroup.Get(_.genegroup, _.avg_read_depth).GroupBy(_.genegroup)())]) / fieldnames
DESEQ_genes = DESEQ_genes[~(_.Get(1).IsMissing() | _.Get(len(DESEQ_genes.Names)-1).IsMissing())].ReplaceMissing(1.0)

Export(DESEQ_genes, '%s/DESEQ_genes.input.tsv' % (output_dir));

###############################################################################
# Groups with chromosomes

Dchromosome = (DESEQ_genes | Match(_.genegroup.Each(lambda x: x.split(',')[0]).Cast(str), _.id) | gffs[0][_.feature == 'mRNA'].Get(_.id, _.seqname)).Without(_.id).GroupBy(_.seqname).Copy()
Dchromosomegroup = Dchromosome
for field in xrange(1,len(Dchromosome.Names)-1):
  print field;
  Dchromosomegroup = Dchromosomegroup.To(field, Do=_.Sum()).Copy()
#efor
DESEQ_chrom = Dchromosomegroup.Get(_.seqname / 'chrgroup', *[_.Get(x) for x in xrange(1, len(Dchromosomegroup.Names)-1)])

Export(DESEQ_chrom, '%s/DESEQ_chrom.input.tsv' % (output_dir));

###############################################################################
# Groups with tissues

Dtissue = DESEQ_chrom.Without(_.chrgroup).Sum()
tissuegroups = [ (x, x+1, y, y+1) for (x,y) in zip(range(0, len(Dtissue.Names)/2,2) , range(len(Dtissue.Names)/2, len(Dtissue.Names), 2))]
Dtissuegroup = []
for group in tissuegroups:
  Dtissuegroup.append((Dtissue.Names[group[0]].split('|')[1],) + Dtissue.Get(*list(group))())
#efor
DESEQ_tissue = Rep(Dtissuegroup[0:]) / (('tissuegroup',) + tuple([ Dtissue.Names[x].split('|')[0] + '|tissue' for x in tissuegroups[0]]))

Export(DESEQ_tissue, '%s/DESEQ_tissue.input.tsv' % (output_dir));

###############################################################################
###############################################################################
###############################################################################
# HERE WE NEED TO RUN DESEQ/EdgeR/whatever for normalization!!!
cmds = []
cmds.append("Rscript deseq.R '%s/DESEQ_genes.input.tsv' %s/DESEQ_genes" % (output_dir, output_dir));
cmds.append("Rscript deseq.R '%s/DESEQ_chrom.input.tsv' %s/DESEQ_chrom" % (output_dir, output_dir));
cmds.append("Rscript deseq.R '%s/DESEQ_tissue.input.tsv' %s/DESEQ_tissue" % (output_dir, output_dir));

utils.run_seq_cmds(cmds)

###############################################################################
###############################################################################
###############################################################################
# Load in DESEQ results


DT_genes  = Read('%s/DESEQ_genes.tests.tsv' % output_dir, delimiter='\t').Detect() / ('condition', 'conda', 'condb', 'id', 'basemean' , 'basemeana', 'basemeanb', 'foldchange', 'log2foldchange', 'pval', 'padj')
DT_chrom  = Read('%s/DESEQ_chrom.tests.tsv' % output_dir, delimiter='\t').Detect() / ('condition', 'conda', 'condb', 'id', 'basemean' , 'basemeana', 'basemeanb', 'foldchange', 'log2foldchange', 'pval', 'padj')
DT_tissue = Read('%s/DESEQ_tissue.tests.tsv' % output_dir, delimiter='\t').Detect() / ('condition', 'conda', 'condb', 'id', 'basemean' , 'basemeana', 'basemeanb', 'foldchange', 'log2foldchange', 'pval', 'padj')

# We have two homokaryons, Homokaryon1 and Homokaryon2, H1 and H2

DT_genes_H1 = (DT_genes.To(_.id, Do=_.Each(lambda x: x.split(',')[0]).Cast(str)) | Match(_.id, jointype='left', merge_same=True) | gffs[0].Get(_.seqname, _.id))
DT_genes_H2 = (DT_genes.To(_.id, Do=_.Each(lambda x: x.split(',')[1]).Cast(str)) | Match(_.id, jointype='left', merge_same=True) | gffs[1].Get(_.seqname, _.id))

all_counted_ids = (DT_genes_H1.id | Stack | DT_genes_H2.id).Copy();

upregulated_in_H1 = DT_genes_H1[_.padj < .05][_.foldchange < 1.0/3.0]
upregulated_in_H2 = DT_genes_H2[_.padj < .05][_.foldchange > 3]

diff_reg_count = (upregulated_in_H1.GroupBy(_.condition).Get(_.condition, _.conda.Count()) | Match(_.condition) | upregulated_in_H2.GroupBy(_.condition).Get(_.condition, _.conda.Count())).Sort(_.condition) / ('condition', 'h1_up', 'h2_up')
Export(diff_reg_count, '%s/diff_reg_count.tsv' % output_dir);

###############################################################################
# Counts of diff ex. genes

diff_regulated = DT_genes[_.padj < .05][(_.foldchange < 1.0/3.0) | (_.foldchange > 3)]

Export(diff_regulated, '%s/diff_regulated.tsv' % output_dir)
Export(diff_regulated.id.Unique(), '%s/diff_regulated_ids.tsv' % output_dir)

# How many?
n_diff_regulated_ids = diff_regulated.id.Unique().Shape()()

# How many up in H1?
diff_regulated[_.foldchange < 1.0/3.0].id.Unique().Shape()()
# How many up in H2?
diff_regulated[_.foldchange > 3.0].id.Unique().Shape()()

# Overlap
set(diff_regulated[_.foldchange < 1.0/3.0].id.Unique()()) & set(diff_regulated[_.foldchange > 3.0].id.Unique()())

# How many in all conditions H1
diff_regulated[_.foldchange < 1.0/3.0].GroupBy(_.id).Get(_.id, _.condition.Count())[_.condition == 17].Shape()

# How many in all conditions H2
diff_regulated[_.foldchange > 3.0].GroupBy(_.id).Get(_.id, _.condition.Count())[_.condition == 17].Shape()

n_in_h1 = upregulated_in_H1.id.Unique().Shape()
n_in_h2 = upregulated_in_H2.id.Unique().Shape()



###############################################################################
###############################################################################
###############################################################################
  # Look for enrichments in these genes

#def format_upregulated_conds(data):
#  conds = data.GroupBy(_.condition).Get(_.condition, _.id, _.seqname)
#  any   = data.To(_.condition, Do=_.Each(lambda x: 'any').Cast(str)).Get(_.condition, _.id, _.seqname).Unique().GroupBy(_.condition)
#  all   = data.GroupBy(_.id).Get(_.id, _.seqname.Count() / 'n', _.seqname[0])[_.n == len(label_names)].GroupBy(_.n).Get(_.n.Each(lambda x: 'all').Cast(str), _.id, _.seqname)
#  return ((conds | Stack | any ) | Stack | all) / ('cond', 'id', 'seqname');
##edef
#
#from scipy import stats as ssp;
#
#def enrich(data):
#  data = data / ('name', 'a', 'b', 'c', 'd')
#  tested = data.Get(_.name, _.a, _.b, _.c, _.d, _.Get(_.a, _.b, _.c, _.d).Each(lambda a,b,c,d: ssp.fisher_exact([ [a, max(b, 1)], [max(c,1), d] ])) / 'test')
#  return tested.Get(_.name, _.a, _.b, _.c, _.d, _.test.Each(lambda x: x[0]).Cast(float) / 'r', _.test.Each(lambda x: x[1]).Cast(float) / 'p').Copy()
##edef
#
#def enriched_chromosome(data):
#  data = data / ('seqname', 'ndiff', 'ngenes')
#  total_diff  = data.ndiff.Sum()()
#  total_genes = data.ngenes.Sum()()
#  enrich_data = data.Get(_.seqname, _.ndiff / 'a', (_.ngenes - _.ndiff) / 'b', _.ndiff.Each(lambda x: total_diff - x).Cast(int) / 'c', _.Get(_.ngenes, _.ndiff).Each(lambda x, y: (total_genes - (x - y))) / 'd')
#  return enrich(enrich_data)
##edef
#
#upregulated_conds_H1 = format_upregulated_conds(upregulated_in_H1)
#upregulated_conds_H2 = format_upregulated_conds(upregulated_in_H2)
#
#all_enrichtests = []
#for org, upregulated_conds, genome in zip(['H1', 'H2'], [ upregulated_conds_H1, upregulated_conds_H2], gffs):
#  for cond in upregulated_conds.cond():
#    print cond
#    data_left  = upregulated_conds.Flat()[_.cond == cond].GroupBy((_.cond, _.seqname)).Get(_.cond[0], _.id.Count() / 'ndiff', _.seqname[0])
#    data_right = genome.GroupBy(_.seqname)[_.feature == 'gene'][_.id.In(all_counted_ids)].Get(_.seqname, _.id.Count() / 'ngenes')
#    data = data_left | Match(_.seqname, jointype='full', merge_same=True) | data_right
#    data = data.Get(_.seqname, _.ndiff, _.ngenes).ReplaceMissing()
#    all_enrichtests.append(enriched_chromosome(data).To(_.name, Do=_.Each(lambda x: org + '|' + cond + '|' + x).Cast(str)).Copy())
#  #efor
##efor
#
#enrichtests = all_enrichtests[0]
#for et in all_enrichtests[1:]:
#  enrichtests = enrichtests | Stack | et;
##efor
#
#Export(enrichtests, '%s/enrichtests.tsv' %output_dir)
#
#enrichtests[_.p < 0.05 / enrichtests.Shape()()].Show()

#scaffold12b_genes = upregulated_conds_H1[_.seqname == 'scaffold_12b'].id.Flat().Unique().Sort().Show()
#mapping0 = Read('%s/mapping_0.fasta' % output_dir)

###############################################################################
###############################################################################
###############################################################################
# Look up PER GENE ratios...

psuedocount = 0.001

sample_pairs = zip(range(1, len(DESEQ_genes.Names)/2 + 1), range(len(DESEQ_genes.Names)/2 + 1,len(DESEQ_genes.Names)+1))
repl_pairs   = zip(range(1, len(sample_pairs) + 1,2), range(2, len(sample_pairs) + 1,2))
pair_names = [ DESEQ_genes.Names[i].split('|')[1] for (i,j) in repl_pairs ]
DI_gene = DESEQ_genes.Get(_.genegroup, *[ ((_.Get(x)+psuedocount) / (_.Get(y)+psuedocount)) for (x,y) in sample_pairs ]).Detect()

  # Calculate average for each replicate condition
combined_groups = [ range(1, len(DESEQ_genes.Names)/2 + 1,2),
                    range(2, len(DESEQ_genes.Names)/2 + 1,2),
                    range(len(DESEQ_genes.Names)/2 + 1, len(DESEQ_genes.Names), 2),
                    range(len(DESEQ_genes.Names)/2 + 2, len(DESEQ_genes.Names), 2) ];

comb_groups_DI = DESEQ_genes.Get(_.genegroup, *[Tuple(*[_.Get(x) for x in y ]).Each(lambda x: np.sum(x)).Cast(int) for y in combined_groups])

  # Calculate average ratio for each replicate
comb_groups_DI = comb_groups_DI.Get(_.genegroup, (_.Get(1)+1.0) / (_.Get(3)+1.0), (_.Get(2)+1.0) / (_.Get(4)+1.0))

  # Add this average to the pairs
repl_pairs.append((max(max(repl_pairs))+1, max(max(repl_pairs))+2))
pair_names.append('all')

  # Attach to DI_gene
DI_gene = DI_gene | Match(_.genegroup, merge_same='equi') | comb_groups_DI

  # Make a huge list of it!
gene_stats = [];
for i, (x, y) in enumerate(repl_pairs):
   print i, x, y
   gene_cond_stats = DI_gene.To(x, Do=_.Each(lambda x: np.log(x)).Cast(float)).To(y, Do=_.Each(lambda x: np.log(x)).Cast(float));
   gene_cond_stats = gene_cond_stats.Get(_.genegroup,
                                        _.genegroup.Each(lambda x: pair_names[i]).Cast(str) / 'cond',
                                        _.Get(x) / 'r1', _.Get(y) / 'r2',
                                       ((_.Get(x) + _.Get(y)) / 2.0) / 'mean')

   gene_cond_stats = gene_cond_stats.Get(_.genegroup, _.cond, _.r1, _.r2, _.mean, ((_.r1 - _.mean) * (_.r1 - _.mean) + (_.r2 - _.mean) * (_.r2 - _.mean) / 2) / 'var')
   gene_cond_stats = gene_cond_stats.Get(_.genegroup, _.cond, _.r1, _.r2, _.mean, _.var.Each(lambda x: np.sqrt(x)).Cast(float) / 'sd', _.var )
   gene_cond_stats = gene_cond_stats.Get(_.genegroup, _.cond, _.r1, _.r2, _.mean, _.var.Each(lambda x: np.sqrt(x)).Cast(float) / 'sd', _.var, (_.sd / np.sqrt(2)) * max(t.interval(0.90, 1)) / 'confidence' ).Copy()
   gene_cond_stats = gene_cond_stats.To(_.mean, Do=_.Each(lambda x: np.exp(x)).Cast(float)).Copy()
   gene_stats.append(gene_cond_stats);
#efor

all_gene_stats = gene_stats[0];
for s in gene_stats[1:]:
  all_gene_stats = all_gene_stats | Stack | s;
#efor
all_gene_stats = all_gene_stats.Copy()

gene_ratio_stats = (all_gene_stats | Match(_.genegroup.Each(lambda x: x.split(',')[0]).Cast(str), _.id,jointype='left') | gffs[0][_.feature == 'mRNA'].Get(_.id, _.seqname, _.start, _.end)).Without(_.id)
gene_ratio_stats = gene_ratio_stats.Copy();

Export(gene_ratio_stats, '%s/gene_ratio_stats.tsv' % output_dir)

genome_chromosomes = genomes[0].Get(_.f0, _.f0.Each(lambda x: 1).Cast(int) /'start', _.seq.Each(lambda x: len(x)).Cast(int) / 'end');
#chrom  chromStart  chromEnd  name    gieStain
ideogram = genome_chromosomes.Get(_.f0, _.start, _.end, _.f0 / 'name', _.f0.Each(lambda x: (1,1,1)).Detect() / 'color');
ideogene = gene_ratio_stats.Get(_.seqname, _.cond, _.start, _.end, _.genegroup, _.mean,
                                _.mean.Each(lambda x: utils.rgb2frac(utils.convert_to_rgb(0.0, 1.0, utils.norm(0.5, 1.0, 2.0, x)))) / 'colors')
ideogene_abs = gene_ratio_stats.Get(_.seqname, _.cond, _.start, _.end, _.genegroup, _.mean,
                                     _.mean.Each(lambda x: utils.rgb2frac(utils.convert_to_rgb(0.0, 1.0, utils.norm(0.9, 1.0, 1.1, x)))) / 'colors')

# Show the genes that we have
gd.draw_genome(utils.natural_sort(genome_chromosomes.f0.Unique()()),
               utils.rep_to_df(ideogram),
               utils.rep_to_df(ideogene.To(_.colors, Do=_.Each(lambda x: (0.0, 0.0, 0.0)))[_.cond == 'all'].Without(_.cond)),
               '%s/genes.differentiable.png' % (output_dir))

# Plot all the genes
for condition in ideogene.cond.Unique()():
  gd2.draw_genome(utils.natural_sort(genome_chromosomes.f0.Unique()()),
                 utils.rep_to_df(ideogram),
                 utils.rep_to_df(ideogene[_.cond == condition].Without(_.cond)),
                 '%s/genes.gene_ratio_chromosome.%s' % (output_dir, condition));
#efor

# Plot all the genes, with a threshold on the color
for condition in ideogene.cond.Unique()():
  gd2.draw_genome(utils.natural_sort(genome_chromosomes.f0.Unique()()),
                 utils.rep_to_df(ideogram),
                 utils.rep_to_df(ideogene_abs[_.cond == condition].Without(_.cond)),
                 '%s/genes.gene_ratio_chromosome.abs.%s' % (output_dir, condition))

#efor


###############################################################################
###############################################################################
###############################################################################
  # DRAW chromosome heatmap
  # Look up AVERAGE CHROMOSOME RATIOS from gene level

def draw_heatmap(allstats, output_prefix='', row_order=None, col_order=None, vmin=None, vmax=None, center=1.0):
  conds, seqname, mean = allstats.Get(_.cond, _.seqname, _.mean)();
  df = pd.DataFrame({ 'conds': conds, 'seqname': seqname, 'mean':mean})
  df = df.pivot('conds', 'seqname', 'mean')
  df = df.reindex(index=row_order) if row_order is not None else df
  df = df.reindex(columns=col_order) if col_order is not None else df
  plt.cla(); plt.clf();
  plt.figure(figsize=(12,6))
  ax = sns.heatmap(df, center=center, vmin=vmin, vmax=vmax)
  plt.savefig('%s.svg' % output_prefix)
  plt.savefig('%s.png' % output_prefix)
#edef

  # Use the gene ratio stats stuff to see which nucleii are more active than the others
DI_chr_ratio = gene_ratio_stats.GroupBy((_.seqname, _.cond)).Get(_.cond[0], _.seqname[0], _.mean.Each(lambda x: np.log(x)).Mean().Each(lambda x: np.exp(x)).Cast(float)).Show()
DI_chr_all   = DI_chr_ratio.GroupBy(_.cond).Get(_.cond, _.seqname.Each(lambda x: 'all').Cast(str)[0], _.mean.Each(lambda x: np.log(x)).Mean().Each(lambda x: np.exp(x)).Cast(float))

DI_chr_ratio_all = DI_chr_ratio | Stack | DI_chr_all

  # Draw without 'all', gives better overview. 'all' displaces mean
row_order = [ x.lower() for x in data_file.Get(_.replicate_id, _.replicate_label).Unique().Sort(_.replicate_id).replicate_label() ]
col_order = ['all'] + utils.natural_sort(DI_chr_ratio.seqname.Unique()())
draw_heatmap(DI_chr_ratio_all, output_prefix = '%s/chr.gene_ratio' % output_dir, row_order=row_order, col_order=col_order)

#Draw another way
DI_chr_gene_logs = gene_ratio_stats.GroupBy((_.seqname, _.cond)).Get(_.cond[0], _.seqname[0], _.mean.Each(lambda x: np.log(x)).Cast(float)).Sort(_.mean)

###############################################################################
  # Look up PER CHROMOSOME READ ratios

DI_chr = DESEQ_genes | Match(_.genegroup.Each(lambda x: x.split(',')[0]).Cast(str), _.id, merge_same='equi', jointype='left') | gffs[0].Get(_.id, _.seqname)
DI_chr = DI_chr.Without(_.id);

DI_chr = DI_chr.GroupBy(_.seqname)

for id, name in enumerate(DI_chr.Names[1:-1], 1):
  DI_chr = DI_chr.To(id, Do=_.Sum())
#efor

DI_chr = DI_chr.Get(_.seqname, *[ x[0] for x in enumerate(DI_chr.Names[1:-1], 1) ]).ReplaceMissing()

conditions = [ x for x in enumerate(DI_chr.Names[1:], 1)]
condition_pairs = zip(conditions[:len(conditions)/2], conditions[len(conditions)/2:])

stats = [];
for ((a_id_rep1, conda_rep1), (b_id_rep1, condb_rep1)), ((a_id_rep2, conda_rep2), (b_id_rep2, condb_rep2)) in zip(condition_pairs[0::2], condition_pairs[1::2]):
  cond = conda_rep1.split('|')[1]
  cond_stats = DI_chr.Get(_.seqname, _.Get(a_id_rep1).Cast(float) / _.Get(b_id_rep1).Cast(float), _.Get(a_id_rep2).Cast(float) / _.Get(b_id_rep2).Cast(float)  )   / ('seqname', 'r1', 'r2');
  cond_stats = cond_stats.To(_.r1, Do=_.Each(lambda x: np.log(x)).Cast(float)).To(_.r2, Do=_.Each(lambda x: np.log(x)).Cast(float));
  cond_stats = cond_stats.Get(_.seqname, _.r1, _.r2, ((_.r1 + _.r2) / 2) / 'mean' )
  cond_stats = cond_stats.Get(_.seqname, _.r1, _.r2, _.mean, ((_.r1 - _.mean) * (_.r1 - _.mean) + (_.r2 - _.mean) * (_.r2 - _.mean) / 2) / 'var')
  cond_stats = cond_stats.Get(_.seqname, _.r1, _.r2, _.mean, _.var.Each(lambda x: np.sqrt(x)).Cast(float) / 'sd', _.var );
  cond_stats = cond_stats.Get(_.seqname, _.r1, _.r2, _.mean, _.sd, _.var, (_.sd / np.sqrt(2)) * max(t.interval(0.90, 1)) / 'confidence' )
  cond_stats = cond_stats.To(_.mean, Do=_.Each(lambda x: np.exp(x)).Cast(float)).Copy()
  stats.append(cond_stats.Get(_.seqname.Each(lambda x: cond).Cast(str) / 'cond', *cond_stats.Names).Copy());
#efor

allstats = stats[0];
for s in stats[1:]:
  allstats = allstats | Stack | s;
#efor

  # Add the data from ALL the chromosomes
condlibratios = Read('%s/DESEQ_genes.condlibratios.tsv' % output_dir).Detect() / ('cond', 'r1', 'r2', 'mean', 'sd', 'var');
condlibratios = condlibratios.Get(_.cond, _.r1, _.r2, _.mean, _.sd, _.var, (_.sd / np.sqrt(2)) * max(t.interval(1 - 0.05/(condlibratios.cond.Shape()()), 1)) / 'confidence')

allstats = condlibratios.Get(_.cond, _.cond.Each(lambda x: 'all').Cast(str) / 'seqname', _.r1, _.r2, _.mean, _.sd, _.var, _.confidence) | Stack | allstats

###############################################################################
  # ACTUALLY DRAW THE HEATMAP

  # Get log values again

log_rr = DI_chr_ratio_all.To(_.mean, Do=_.Each(lambda x: np.log2(x)).Cast(float))
log_gr = allstats.To(_.mean, Do=_.Each(lambda x: np.log2(x)).Cast(float))

  # DETERMINE VMIN VMAX

minval_rr, maxval_gr = DI_chr_ratio_all.Get(_.mean.Min(), _.mean.Max())()
minval_gr, maxval_rr = allstats.Get(_.mean.Min(), _.mean.Max())()

minmin = min(minval_rr, minval_gr)
maxmax = max(maxval_gr, maxval_rr);

maxval = max(maxmax, 1.0/minmin)
minval = min(minmin, 1.0/maxmax)

minlogval_rr, maxlogval_rr = log_gr.Get(_.mean.Min(), _.mean.Max())()
minlogval_gr, maxlogval_gr = log_rr.Get(_.mean.Min(), _.mean.Max())()

minlogval = min(minlogval_rr, minlogval_gr, -maxlogval_rr, -maxlogval_gr)
maxlogval = -minlogval


  # GENE RATIO
row_order = [ x.lower() for x in data_file.Get(_.replicate_id, _.replicate_label).Unique().Sort(_.replicate_id).replicate_label() ]
col_order = ['all'] + utils.natural_sort(DI_chr_ratio.seqname.Unique()())
draw_heatmap(DI_chr_ratio_all, output_prefix = '%s/chr.gene_ratio' % output_dir, row_order=row_order, col_order=col_order, vmin=minval,  vmax=maxval)


  # READ RATIO
row_order = [ x.lower() for x in data_file.Get(_.replicate_id, _.replicate_label).Unique().Sort(_.replicate_id).replicate_label() ]
col_order = utils.natural_sort(allstats.seqname.Unique()())
draw_heatmap(allstats, '%s/chr.read_ratio' % output_dir, row_order=row_order, col_order=col_order, vmin=minval,  vmax=maxval)



  # GENE RATIO LOG
row_order = [ x.lower() for x in data_file.Get(_.replicate_id, _.replicate_label).Unique().Sort(_.replicate_id).replicate_label() ]
col_order = ['all'] + utils.natural_sort(DI_chr_ratio.seqname.Unique()())
draw_heatmap(log_rr,'%s/chr.gene_ratio_log' % output_dir, row_order=row_order, col_order=col_order, vmin=minlogval, vmax=maxlogval, center=0.0)
  
  
  # READ RATIO LOG
row_order = [ x.lower() for x in data_file.Get(_.replicate_id, _.replicate_label).Unique().Sort(_.replicate_id).replicate_label() ]
col_order = utils.natural_sort(allstats.seqname.Unique()())
draw_heatmap(log_gr, '%s/chr.read_ratio_log' % output_dir, row_order=row_order, col_order=col_order, vmin=minlogval, vmax=maxlogval, center=0.0)

  # One sample t-test to check for dominance of one chromosome or the other
import scipy
scipy.stats.ttest_1samp(log_gr.mean(), 0)

###############################################################################
###############################################################################
###############################################################################
  #Per tissue ratios

# Gene level ratio
DI_tissue_ratio = DI_chr_ratio | Stack | DI_chr_all
                                                                                                                                                                                                                                             # Read level ratio
allstats = allstats

# Determine maximum values
minval_gr, maxval_gr = DI_tissue_ratio.Get(_.mean.Min(), _.mean.Max())()
minval_rr, maxval_rr = allstats.Get(_.mean.Min(), _.mean.Max())()

minmin = min(minval_rr, minval_gr)
maxmax = max(maxval_gr, maxval_rr);

maxval = max(maxmax, 1.0/minmin)
minval = min(minmin, 1.0/maxmax)

DI_tissue_ratio_log = DI_tissue_ratio.To(_.mean, Do=_.Each(lambda x: np.log2(x)).Cast(float))
allstats_log        = allstats.To(_.mean, Do=_.Each(lambda x: np.log2(x)).Cast(float))

minval_grl, maxval_grl = DI_tissue_ratio_log.Get(_.mean.Min(), _.mean.Max())()
minval_rrl, maxval_rrl = allstats_log.Get(_.mean.Min(), _.mean.Max())()

minmin_log = min(minval_grl, minval_rrl)
maxmax_log = max(maxval_grl, maxval_rrl)

maxval_log = max(maxmax_log, -minmin_log)
minval_log = -maxval_log

###############################################################################



  # Look up average tissue ratios from gene level
DI_tissue_ratio_norm = DI_tissue_ratio.To(_.mean, Do=_.Each(lambda x: utils.norm(minval, 1.0, maxval, x)).Cast(float))

def minmaxnorm(min, max, x):
  return float(x - min) / float(max - min)
#edef

cmap = matplotlib.cm.get_cmap('RdBu_r')

DI_tissue_ratio_norm_log = DI_tissue_ratio_log.To(_.mean, Do=_.Each(lambda x: minmaxnorm(minval_log, maxval_log, x)).Cast(float))

# Draw GENE RATIOS
colors = [(0, 0, 255), (255, 255, 255), (255, 0, 0)]  # [BLUE, WHITE, RED]
#colors = [ (5,48,97), (255,255,255), (103,30,31)] # [NAVYBLUE,WHITE,WINERED]
defaults = { '!!H1!!': 'P1', '!!H2!!': 'P2', '!!ratio_min!!': '%1.1f' % minval, '!!ratio_max!!': '%1.1f' % maxval, '!!ratio_mid!!' : '1' }
if tissue_figure is not None:
   for chrid in DI_tissue_ratio_norm.seqname.Unique()():
     replacevals = defaults
     replacevals.update(dict(zip(*DI_tissue_ratio_norm[_.seqname == chrid].Get(_.cond.Each(lambda x: '!!' + x + '!!').Cast(str),
                                                                              _.mean.Each(lambda x: utils.rgb2hex(utils.convert_to_rgb(0, 1, x, colors))))())))
     svg = svg_utils.read_svg(tissue_figure);
     svg = svg_utils.replacem(svg, replacevals);
     svg_utils.write_svg(svg, '%s/tissue.gene_ratio.%s.svg' % (output_dir, chrid))
   #efor
#fi

# Draw GENE RATIOS (LOG)
colors = [(0, 0, 255), (255, 255, 255), (255, 0, 0)]  # [BLUE, WHITE, RED]
#colors = [ (5,48,97), (255,255,255), (103,30,31)] # [NAVYBLUE,WHITE,WINERED]
defaults = { '!!H1!!': 'P1', '!!H2!!': 'P2', '!!ratio_min!!': '%1.1f' % minval_log, '!!ratio_max!!': '%1.1f' % maxval_log, '!!ratio_mid!!' : '0' }
if tissue_figure is not None:
   for chrid in DI_tissue_ratio_norm.seqname.Unique()():
     replacevals = defaults
     replacevals.update(dict(zip(*DI_tissue_ratio_norm[_.seqname == chrid].Get(_.cond.Each(lambda x: '!!' + x + '!!').Cast(str), _.mean.Each(lambda x: matplotlib.colors.rgb2hex(cmap(x))[1:]).Cast(str))())))
     svg = svg_utils.read_svg(tissue_figure);
     svg = svg_utils.replacem(svg, replacevals);
     svg_utils.write_svg(svg, '%s/tissue.gene_ratio_log.%s.svg' % (output_dir, chrid))
   #efor
#fi

###############################################################################
  # Draw READ ratios

colors = [(0, 0, 255), (255, 255, 255), (255, 0, 0)]  # [BLUE, WHITE, RED]
#colors = [ (5,48,97), (255,255,255), (103,30,31)] # [NAVYBLUE,WHITE,WINERED]
defaults = { '!!H1!!': 'P1', '!!H2!!': 'P2', '!!ratio_min!!': '%1.1f' % minval, '!!ratio_max!!': '%1.1f' % maxval, '!!ratio_mid!!' : '1' }
allstats_mean_norm     = allstats.Get(_.mean.Each(lambda x: utils.norm(minval, 1.0, maxval, x)).Cast(float) / 'norm', *allstats.Names)
allstats_mean_norm_log = allstats_log.To(_.mean, Do=_.Each(lambda x: minmaxnorm(minval_log, maxval_log, x)).Cast(float))

if tissue_figure is not None:
   reload(svg_utils)
   reload(utils);
   for chrid in allstats.seqname.Unique()():
     replacevals = defaults
     replacevals.update(dict(zip(*allstats_mean_norm[_.seqname == chrid].Get(_.cond.Each(lambda x: '!!' + x + '!!').Cast(str),
                                                                            _.norm.Each(lambda x: utils.rgb2hex(utils.convert_to_rgb(0, 1, x, colors))))())))
     svg = svg_utils.read_svg(tissue_figure);
     svg = svg_utils.replacem(svg, replacevals);
     svg_utils.write_svg(svg, '%s/tissue.read_ratio.%s.svg' % (output_dir, chrid))
   #efor
#fi

defaults = { '!!H1!!': 'P1', '!!H2!!': 'P2', '!!ratio_min!!': '%1.1f' % minval_log, '!!ratio_max!!': '%1.1f' % maxval_log, '!!ratio_mid!!' : '0' }
if tissue_figure is not None:
   reload(svg_utils)
   reload(utils);
   for chrid in allstats.seqname.Unique()():
     replacevals = defaults
     replacevals.update(dict(zip(*allstats_mean_norm_log[_.seqname == chrid].Get(_.cond.Each(lambda x: '!!' + x + '!!').Cast(str), _.mean.Each(lambda x: matplotlib.colors.rgb2hex(cmap(x))[1:]).Cast(str))())))
     svg = svg_utils.read_svg(tissue_figure);
     svg = svg_utils.replacem(svg, replacevals);
     svg_utils.write_svg(svg, '%s/tissue.read_ratio_log.%s.svg' % (output_dir, chrid))
   #efor
#fi


###############################################################################
###############################################################################
###############################################################################
  # CLUSTER ANALYSIS?
#def wisecondor(signal,  nbins=10, thresh=3.0):
#
#  mu = np.mean(signal);
#  sd = np.std(signal);
#
#  z_signal = (signal - mu) / sd
#
#  sl_window = np.array([ np.sum(z_signal[i-nbins:i+nbins+1]) / (nbins*2 +1) for i in xrange(nbins, len(signal) - nbins) ])
#
#  return sl_window, np.abs(sl_window) >= thresh, mu, sd;
##edef
#
#def wisecondor_called_regions(sl_window, sl_thresh, nbins=10):
#  in_sig_region = False
#  sig_region_start = 0;
#  sig_regions = []
#  for i, sig in enumerate(sl_thresh):
#    if sig and not(in_sig_region):
#      sig_region_start = i + nbins;
#    elif not(sig) and in_sig_region:
#      sig_regions.append((np.mean(sl_window[sig_region_start-nbins:i+1]), sig_region_start, i-1+2*nbins));
#    #fi
#    in_sig_region = sig;
#  #efor
#  return sig_regions;
##edef
#
#nbins=5;
#thresh=2;
#found_regions = []
#for condition in gene_ratio_stats.cond.Unique()():
#  for scaffold in gene_ratio_stats.seqname.Unique()():
#    chr_genes, chr_ratios = gene_ratio_stats[_.cond == condition][_.seqname == scaffold].Sort(_.start).Get(_.genegroup, _.mean)();
#    wisecondor_sl_window, wisecondor_thresh, wisecondor_mu, wisecondor_std = wisecondor(chr_ratios, nbins, thresh)
#    [ found_regions.append((condition, scaffold) + tuple(x))for x in wisecondor_called_regions(wisecondor_sl_window, wisecondor_thresh, nbins) ]
#  #efor
##efor
#
def simple_grouping(signal, min_cluster_length=10):
  def sign(v):
    return v >= 1.0
  #edef

  clusters = [];

  last_sign = sign(signal[0]);

  in_cluster = False
  start_cluster = 0;
  for i, v in enumerate(signal,1):
    sv = sign(v)
    if not(in_cluster) and (sv == last_sign):
      in_cluster = True;
      start_cluster = i-1;
    elif in_cluster and sv != last_sign:
      in_cluster = False;
      if i-1-start_cluster > min_cluster_length:
        clusters.append((start_cluster, i-1, last_sign));
      #fi
    #fi
    last_sign = sv;
  #efor

  return clusters;
#edef

found_clusters = []
for condition in gene_ratio_stats.cond.Unique()():
  for scaffold in gene_ratio_stats.seqname.Unique()():
    chr_genes, chr_ratios, chr_starts, chr_ends = gene_ratio_stats[_.cond == condition][_.seqname == scaffold].Sort(_.start).Get(_.genegroup, _.mean, _.start, _.end)();
    [ found_clusters.append((condition, scaffold, chr_starts[s], chr_ends[e], sign)) for (s, e, sign) in simple_grouping(chr_ratios) ]
  #efor
#efor

###############################################################################
###############################################################################
###############################################################################

# Question: So now we know that P1 is more active on average at a GENE level, why is it that P2 is more active at a READ level?
# Let's look at the genes which are very HIGHly active, in general(where the GENE level ratio isn't affected so much, but the READ level is...


DT_gene_chr = (DT_genes | Match(_.id.Each(lambda x: x.split(',')[0]).Cast(str), _.id, merge_same='equi', jointype='left') | gffs[0].Get(_.id, _.seqname)).Without(_.id.R)

  # Cond A Greater Than B
#AGTB = DT_gene_chr.GroupBy((_.condition, _.seqname))[_.basemeanb < _.basemeana][_.padj > 0.05].Get(_.condition[0], _.seqname[0], _.id.Count() / 'agtb_count')
  # Cond B Greated Than A
#BGTA = DT_gene_chr.GroupBy((_.condition, _.seqname))[_.basemeana < _.basemeanb][_.padj > 0.05].Get(_.condition[0], _.seqname[0], _.id.Count() / 'bgta_count')

#(AGTB | Match((_.condition, _.seqname), merge_same='equi') | BGTA).Without(_.condition.R, _.seqname.R)

def determine_read_ratio_plot(p1_expr, p2_expr):

  psuedocount = 0.0001
  maxvals = np.array([ max(a,b) for (a,b) in zip(p1_expr, p2_expr)]);
  p1_sort, p2_sort, maxval_sort = zip(*sorted(zip(p1_expr, p2_expr, maxvals), key=lambda x: x[2]))
  p1_sort = np.array(p1_sort) + psuedocount;
  p2_sort = np.array(p2_sort) + psuedocount;

  p1_sum = psuedocount;
  p2_sum = psuedocount;
  ratio_arr = [];

  for p1, p2 in zip(p1_sort, p2_sort):
    p1_sum += p1;
    p2_sum += p2;

    ratio_arr.append(p1_sum / p2_sum);
  #efor

  return np.array(ratio_arr);
#edef

def determine_gene_ratio_plot(p1_expr, p2_expr):

  psuedocount = 0.0001
  maxvals = np.array([ max(a,b) for (a,b) in zip(p1_expr, p2_expr)]);
  p1_sort, p2_sort, maxval_sort = zip(*sorted(zip(p1_expr, p2_expr, maxvals), key=lambda x: x[2]))
  p1_sort = np.array(p1_sort) + psuedocount;
  p2_sort = np.array(p2_sort) + psuedocount;

  ratio_arr = [];

  for i in xrange(len(p1_sort)):
    ratio_arr.append(np.exp(np.mean([np.log(x) - np.log(y) for (x,y) in zip(p1_sort[:i+1], p2_sort[:i+1])])))
  #efor

  return ratio_arr;

#edef

def draw_scatterplots(D, condition, outprefix):
  plt.cla(); plt.clf();
  f, axes = plt.subplots(4,4, figsize=(40,40));
  axes = [ ax for axl in axes for ax in axl ];
  plt.title(condition);

  Ds = D[_.condition == condition].GroupBy(_.seqname).Get(_.seqname, _.basemeana, _.basemeanb);

  for i, (seqname, p1_expr, p2_expr) in enumerate(zip(*Ds.Sort(_.seqname.Each(lambda x: x.split('_')[-1]).Cast(int))())):
    #axes[i].set(xscale="log", yscale="log")
    sns.regplot(x=p1_expr, y=p2_expr, ci=False, ax=axes[i], fit_reg=False);
    axes[i].set_title('%s (%d) %d above; %d below; %0.3f' % (seqname, len(p1_expr), np.sum(p1_expr > p2_expr), np.sum(p1_expr <= p2_expr), np.sum(p1_expr) / np.sum(p2_expr)));
    axes[i].set_xlabel('P1');
    axes[i].set_ylabel('P2');
    maxmax = max( max(p1_expr), max(p2_expr)) + 500
    axes[i].plot([0,maxmax], [0, maxmax]);
    otherx = axes[i].twinx();
    otherx.plot(sorted(p1_expr), determine_read_ratio_plot(p1_expr, p2_expr), c='r');
    otherx.plot(sorted(p1_expr), determine_gene_ratio_plot(p1_expr, p2_expr), c='g');
    otherx.set_ylim([0,3])
    axes[i].set_xlim([0, maxmax]);
    axes[i].set_ylim([0, maxmax]);
  #efor

  plt.savefig('%s.png' % outprefix);
  plt.savefig('%s.svg' % outprefix);

#edef

for condition in DT_gene_chr.condition.Unique()():
  draw_scatterplots(DT_gene_chr, condition, '%s/dist.read_dist.%s' % (output_dir, condition));
#efor

#OK, so which genes are they??
sorted_extreme_genes = DT_gene_chr.GroupBy(_.id).Get(_.id, _.seqname[0], _.condition, (_.basemeanb - _.basemeana).Cast(int) / 'diff')[_.diff == _.diff.Max()][_.diff.Count() == 1].Flat().GroupBy(_.seqname).Sort(_.diff, descend=True).Sort(_.diff.Mean())

extreme_contributions = DT_gene_chr.GroupBy((_.condition, _.seqname)).Get(_.condition[0], _.seqname[0], _.id, _.basemeana.Cast(float) / (_.basemeana.Sum() + _.basemeanb.Sum()) / 'h1', _.basemeanb.Cast(float) / (_.basemeana.Sum() + _.basemeanb.Sum()) / 'h2')

most_extreme_contributions = extreme_contributions.Sort(_.h2, descend=True).Sort(_.h2[0], descend=True)[(_.h2 > 0.1) | (_.h1 > 0.1)][_.id.Count() > 0].Show()

most_extreme_contributions_P1 = most_extreme_contributions[_.h1 > _.h2].Flat()
most_extreme_contributions_P2 = most_extreme_contributions[_.h2 > _.h1].Flat()

Export(most_extreme_contributions.Flat().Sort(_.seqname, _.id,_.h2 - _.h1, descend=True).Show(), '%s/diff.most_extreme.tsv' % output_dir)


most_extreme_tissue  = Read('A15/diff.most_extreme.tsv').Detect()
most_extreme_compost = Read('A15_compost_nodup/diff.most_extreme.tsv').Detect()

(most_extreme_tissue | Stack | most_extreme_compost).id.Unique()

###############################################################################
###############################################################################
###############################################################################
  # Functional analysis

annotation_file = 'agabi_annotations.tsv';
A  = Read(annotation_file) / ('annot', 'p1_file', 'p2_file');

def cond_func(AD, diff_regulated, condition, parent=None):


  DRC = diff_regulated[_.condition == condition]
  AD  = [ (a,
           b[_.id.In(diff_regulated.id.Each(lambda x: x.split(',')[0]).Cast(str))].To(_.id, Do=_.TakeFrom(diff_regulated.Get(_.id.Each(lambda x: x.split(',')[0]).Cast(str), _.id))),
           c[_.id.In(diff_regulated.id.Each(lambda x: x.split(',')[1]).Cast(str))].To(_.id, Do=_.TakeFrom(diff_regulated.Get(_.id.Each(lambda x: x.split(',')[1]).Cast(str), _.id))))
           for (a,b,c) in AD ]


  upregulated_in_H1 = DRC[_.foldchange < 1.0/3.0].Get(_.id);
  upregulated_in_H2 = DRC[_.foldchange > 3.0].Get(_.id);

  if parent is not None:
    DRP = diff_regulated[_.condition == parent];
    parent_upregulated_in_H1 = DRP[_.foldchange < 1.0/3.0].Get(_.id);
    parent_upregulated_in_H2 = DRP[_.foldchange > 3.0].Get(_.id);

    no_switch_H1 = upregulated_in_H1[_.id.In(parent_upregulated_in_H1)];
    no_switch_H2 = upregulated_in_H2[_.id.In(parent_upregulated_in_H2)];
  #fi

  combined_funcs = AD[0];
  for (a,b,c) in AD[1:]:
    combined_funcs = ('other', (combined_funcs[1] | Stack | b).Unique(), (combined_funcs[2] | Stack | c).Unique())
  #efor

  ADIFF = [ (condition, a,
             b[_.id.In(upregulated_in_H1)].Shape()(), # Number that are annotated as this, and UPREGULATED in H1
             c[_.id.In(upregulated_in_H2)].Shape()(), # Number that are annotated as this, and UPREGULATED in H2
             b[_.id.In(no_switch_H1)].Shape()() if parent is not None else 0,
             c[_.id.In(no_switch_H2)].Shape()() if parent is not None else 0,
            ) for (a,b,c) in AD]
  ADIFF.append( (condition, 'other',
                 upregulated_in_H1.Shape()() - combined_funcs[1][_.id.In(upregulated_in_H1)].Shape()(),
                 upregulated_in_H2.Shape()() - combined_funcs[2][_.id.In(upregulated_in_H2)].Shape()(),
                 no_switch_H1.Shape()()      - sum([ a[4] for a in ADIFF]) if parent is not None else 0,
                 no_switch_H2.Shape()()      - sum([ a[5] for a in ADIFF]) if parent is not None else 0 ));
  ADIFF.append( (condition, 'any',
                 upregulated_in_H1.Shape()(),
                 upregulated_in_H2.Shape()(),
                 no_switch_H1.Shape()() if parent is not None else 0,
                 no_switch_H2.Shape()() if parent is not None else 0));
  return ADIFF

#edef

def traverse_func_tree(AD, diff_regulated, tree, root):

  stack = [ (root, None) ];
  cond_func_table = [];

  while len(stack) > 0:
    node, parent = stack.pop();
    stack.extend([ (c, node) for c in tree[node]] if node in tree else []);
    cond_func_table.extend(cond_func(AD, diff_regulated, node, parent));
  #ewhile

  cond_func_rep = Rep(cond_func_table) / ('condition', 'func', 'h1_up', 'h2_up', 'h1_parent_overlap', 'h2_parent_overlap');

  return cond_func_rep;

#edef


AD = [ (annot, Read(p1_file, delimiter=',')[_.Get(0).In(M.Get(0))] / ('id',), Read(p2_file, delimiter=',')[_.Get(0).In(M.Get(1))] / ('id',)) for (annot, p1_file, p2_file) in zip(*A())]

differentiation = { 'veg':              ['initials'],
                    'initials':         ['ps_stipe_center', 'ps_stipe_shell'],
                    'ps_stipe_center':  ['dif_stipe_center', 'dif_stipe_shell', 'dif_stipe_skin'],
                    'ps_stipe_shell':   ['dif_cap_skin', 'dif_cap_tissue', 'dif_gill_tissue'],
                    'dif_stipe_center': ['yfb_stipe_center'],
                    'dif_stipe_shell':  ['yfb_stipe_shell'],
                    'dif_stipe_skin':   ['yfb_stipe_skin'],
                    'dif_cap_skin':     ['yfb_cap_skin'],
                    'dif_cap_tissue':   ['yfb_cap_tissue'],
                    'dif_gill_tissue':  ['yfb_gill_tissue', 'yfb_veil'] }

compost_diff =  { 'day16': ['pinning'],
                  'pinning': ['flush1'],
                  'flush1': ['postflush1'],
                  'postflush1': ['flush2'],
                  'flush2': ['postflush2'] };


sample_tree, sample_root = (differentiation, 'veg') if tissue_figure is not None else (compost_diff, 'day16')

conf_func_rep = traverse_func_tree(AD, diff_regulated, sample_tree, sample_root)

# Calculate for all samples
cf_all = Rep(cond_func(AD, diff_regulated.To(_.condition, Do=_.Each(lambda x: 'all').Cast(str)).Unique(_.id), 'all', None)) / ('condition', 'func', 'h1_up', 'h2_up', 'h1_parent_overlap', 'h2_parent_overlap')

###############################################################################
###############################################################################
###############################################################################

# Named genes
named = Read('/tudelft.net/staff-groups/ewi/insy/DBL/thiesgehrmann/w/phd/data/A15_homokaryons/names/named_genes_mapping.tsv', delimiter='\t') / ('p1', 'p2', 'agabi2', 'name', 'evalue', 'pident');

upregulated_in_H1[_.id.In(named.p1)]
upregulated_in_H2[_.id.In(named.p2)]

# Draw figure for MnP

mnp_expr         = gene_ratio_stats[_.genegroup == 'AgabiA15p1|4008,AgabiA15p2|3983'].Show().Get(_.cond, _.mean.Each(lambda x: np.log2(x)).Cast(float));
mnp_min, mnp_max = mnp_expr.Get(_.mean.Min(), _.mean.Max())();
min_mnp_log = min(mnp_min, -mnp_max)
max_mnp_log = -min_mnp_log
mnp_expr_norm    = mnp_expr.Get(_.cond,_.mean.Each(lambda x: minmaxnorm(min_mnp_log, max_mnp_log, x)).Cast(float) / 'norm')

defaults         = { '!!H1!!': 'P1', '!!H2!!': 'P2', '!!ratio_min!!': '%1.1f' % min_mnp_log, '!!ratio_max!!': '%1.1f' % max_mnp_log, '!!ratio_mid!!' : '0' }
if tissue_figure is not None:
  reload(svg_utils)
  reload(utils);
  replacevals = defaults
  replacevals.update(dict(zip(*mnp_expr_norm.Get(_.cond.Each(lambda x: '!!' + x + '!!').Cast(str), _.norm.Each(lambda x: matplotlib.colors.rgb2hex(cmap(x))[1:]).Cast(str))())))
  svg = svg_utils.read_svg(tissue_figure);
  svg = svg_utils.replacem(svg, replacevals);
  svg_utils.write_svg(svg, '%s/indiv_genes.mnp_log.svg' % (output_dir));
#fi

mnp_raw_expr = DT_genes[_.id == 'AgabiA15p1|4008,AgabiA15p2|3983'].Get(_.condition, _.basemeana, _.basemeanb, (_.basemeana + _.basemeanb) / 'total')

plt.cla(); plt.clf();
plt.figure(figsize=(12,6))
plt.plot(range(len(mnp_raw_expr.basemeana())), mnp_raw_expr.basemeana(), c='r', label='P1');
plt.plot(range(len(mnp_raw_expr.basemeana())), mnp_raw_expr.basemeanb(), c='b', label='P2');
plt.plot(range(len(mnp_raw_expr.basemeana())), mnp_raw_expr.total(), c='g', label='Total');

plt.xticks(xrange(6), mnp_raw_expr.condition(), rotation='vertical')

plt.savefig('%s/indiv_genes.mnp.expr.png' % output_dir)
plt.savefig('%s/indiv_genes.mnp.expr.svg' % output_dir)



###############################################################################
###############################################################################
###############################################################################
  #PCA plot of different samples
from sklearn.decomposition import PCA

norm_expr       = ((DT_genes.Get(_.conda, _.basemeana.Cast(float)) | Stack | DT_genes.Get(_.condb, _.basemeanb.Cast(float))) / ('condition', 'expr')).GroupBy(_.condition)
norm_expr_names = norm_expr.condition();
norm_expr       = norm_expr.expr()
norm_expr       = np.array([ x for x in norm_expr])

pca = PCA(n_components=3);
DR  = pca.fit_transform(norm_expr)

def plot_pca(pca, DR, names, component1, component2, output_dir):

  plt.cla(); plt.clf();
  plt.figure(figsize=(12,6))
  X = [0] * DR.shape[0] if component1 is None else DR[:,component1]
  Y = [0] * DR.shape[0] if component2 is None else DR[:,component2]
  C = [ 'red' if 'p1' in x else 'blue' for x in names]

  var_expl_1 = 0 if component1 is None else pca.explained_variance_ratio_[component1]
  var_expl_2 = 0 if component2 is None else pca.explained_variance_ratio_[component2]

  component1 = '%d' % (component1+1) if component1 is not None else 'None'
  component2 = '%d' % (component2+1) if component2 is not None else 'None'

  plt.scatter(X,Y, c=C)
  plt.xlabel('Component %s (%f)' % (component1, var_expl_1));
  plt.ylabel('Component %s (%f)' % (component2, var_expl_2));
  fontdict = { 'size': 5, 'ha': 'left', 'va': 'bottom'}
  [ plt.text(x,y, '%s' % norm_expr_names[i], fontdict=fontdict, rotation=90 if (component1 == 'None' or component2 == 'None') else 45) for i,(x,y) in enumerate(zip(X,Y))]

  plt.savefig('%s/pca.samples.%s,%s.png' % (output_dir, component1, component2));
  plt.savefig('%s/pca.samples.%s,%s.svg' % (output_dir, component1, component2));

#edef

plot_pca(pca, DR, norm_expr_names, 0, None, output_dir);
plot_pca(pca, DR, norm_expr_names, 1, None, output_dir);
plot_pca(pca, DR, norm_expr_names, 2, None, output_dir);
plot_pca(pca, DR, norm_expr_names, 0, 1, output_dir);
plot_pca(pca, DR, norm_expr_names, 1, 2, output_dir);

###############################################################################
###############################################################################
###############################################################################
  # Differential protein annotation usage

H1_pfam = Read('/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/annotations/AgabiA15p1.PFAM.txt', delimiter='\t');
H2_pfam = Read('/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/annotations/AgabiA15p2.PFAM.txt', delimiter='\t');

H1_pfam

diff_A = (M | Match(_.f0, _.proteinid, jointype='left') | H1_pfam.GroupBy(_.proteinid).Get(_.proteinid, _.pfam.Unique() / 'h1_pfam')).Without(_.proteinid)
diff_A = (diff_A | Match(_.f1, _.proteinid, jointype='left') | H2_pfam.GroupBy(_.proteinid).Get(_.proteinid, _.pfam.Unique() / 'h2_pfam')).Without(_.proteinid)

domain_features = diff_A.Get(_.f0, _.f1,
                             _.Get(_.h1_pfam, _.h2_pfam).Array().Each(lambda x,y: len(set(x) & set(y))).Cast(int) / 'intersection',
                             _.Get(_.h1_pfam, _.h2_pfam).Array().Each(lambda x,y: len(set(x) | set(y))).Cast(int) / 'union',
                             _.Get(_.h1_pfam, _.h2_pfam).Array().Each(lambda x,y: min(len(x), len(y))).Cast(int) / 'min');

intersection_union = domain_features.GroupBy((_.intersection, _.union)).Get(_.intersection[0], _.union[0], _.f0.Count())

df_domain_features = utils.rep_to_df(intersection_union);
df_domain_features['f0'] = df_domain_features['f0'].astype('int')
df_domain_features = df_domain_features.pivot('intersection', 'union', 'f0')
df_domain_features = df_domain_features.fillna(0)

plt.cla(); plt.clf();
rc={'axes.labelsize': 20, 'font.size': 10, 'legend.fontsize': 15, 'axes.titlesize': 20, 'xtick.labelsize': 15, 'ytick.labelsize': 15}
sns.set_style("whitegrid");
sns.set(rc=rc);
sns.heatmap(df_domain_features[-1::-1], annot=True, fmt='0.0f', cmap="Greys", cbar=False)
plt.savefig('%s/func.union,intersection.png' % output_dir, transparent=True);
plt.savefig('%s/func.union,intersection.svg' % output_dir, transparent=True);

###############################################################################
###############################################################################
###############################################################################
  # Overlap of differential HSE compost & tissue

diff_reg_compost = Read('A15_compost_nodup/diff_regulated.tsv', delimiter='\t').Cast(str, str, str, str, float, float, float, float, float, float, float)
diff_reg_tissue  = Read('A15/diff_regulated.tsv', delimiter='\t').Cast(str, str, str, str, float, float, float, float, float, float, float)

diff_reg_overlap_ids = diff_reg_compost.id.Unique()[_.id.In(diff_reg_tissue.id.Unique())]

diff_reg_overlap = (diff_reg_compost | Stack | diff_reg_tissue)[_.id.In(diff_reg_overlap_ids)]

diff_reg_overlap_H1 = diff_reg_overlap[_.foldchange < 1.0/3.0]
diff_reg_overlap_H2 = diff_reg_overlap[_.foldchange > 3]

diff_reg_overlap[_.id.In(diff_reg_overlap_H2.id)]

  # Which are differentially regulated between compost and tissue (i.e. switching up in tissue or down in tissue, when it was always up in compost)
diff_reg_diff = diff_reg_overlap.GroupBy(_.id)[_.foldchange.Array().Each(lambda x: np.any(x > 3.0) & np.any(x < 1.0/3.0)).Cast(bool)].Get(_.id, _.condition, _.basemeana, _.basemeanb)
diff_reg_diff_table = diff_reg_overlap.GroupBy(_.id)[_.foldchange.Array().Each(lambda x: np.any(x > 3.0) & np.any(x < 1.0/3.0)).Cast(bool)].Get(_.id, _.condition, _.basemeana, _.basemeanb, _.log2foldchange).Flat();

Export(diff_reg_diff_table, '%s/diff.overlap_diff.tsv' % output_dir)


OV_ADIFF = [ (a,
             b[_.id.In(diff_reg_overlap_H1.id.Each(lambda x: x.split(',')[0]).Cast(str))].Shape()(), # Number that are annotated as this, and UPREGULATED in H1
             c[_.id.In(diff_reg_overlap_H2.id.Each(lambda x: x.split(',')[1]).Cast(str))].Shape()(), # Number that are annotated as this, and UPREGULATED in H2
            ) for (a,b,c) in AD]

OV_ADIFF = [ (a,
             b[_.id.In(diff_reg_overlap_H1.id.Each(lambda x: x.split(',')[0]).Cast(str))](), # Number that are annotated as this, and UPREGULATED in H1
             c[_.id.In(diff_reg_overlap_H2.id.Each(lambda x: x.split(',')[1]).Cast(str))](), # Number that are annotated as this, and UPREGULATED in H2
            ) for (a,b,c) in AD]

# Differential domain annotations for these overlapping genes
domain_features.Get(_.f0 + ',' + _.f1, _.intersection, _.union)[_.result.In(diff_reg_overlap_ids)].GroupBy((_.intersection, _.union)).Get(_.intersection[0], _.union[0], _.result.Count())[_.union != _.intersection]

###############################################################################
###############################################################################
###############################################################################
  # KEGG ANNOTS!!!

import imp

load_kegg = imp.load_source('module.name', '/home/nfs/thiesgehrmann/groups/w/phd/tasks/kegg_predictions/load_kegg.py')

diff_reg_all = diff_reg_compost | Stack | diff_reg_tissue

diff_reg_all_H1 = diff_reg_all[_.foldchange < 1.0/3.0]
diff_reg_all_H2 = diff_reg_all[_.foldchange > 3.0]

KEGG_P1 = load_kegg.read_kegg('/home/nfs/thiesgehrmann/groups/w/phd/tasks/kegg_predictions/P1_hier_MODULE.keg') | Match(_.e_query_id, _.f0, merge_same=True) | M
KEGG_P1 = KEGG_P1.Get(_.d, _.d_desc, _.d_path, _.e_ec, _.e_desc, (_.e_query_id + ',' + _.f1) / 'genegroup');
KEGG_P2 = load_kegg.read_kegg('/home/nfs/thiesgehrmann/groups/w/phd/tasks/kegg_predictions/P2_hier_MODULE.keg') | Match(_.e_query_id, _.f1, merge_same=True) | M
KEGG_P2 = KEGG_P2.Get(_.d, _.d_desc, _.d_path, _.e_ec, _.e_desc, (_.f0 + ',' + _.e_query_id) / 'genegroup');

KEGG = KEGG_P1 | Stack | KEGG_P2

KEGG_pathways_P1 = KEGG[_.genegroup.In(diff_reg_all_H1.id)].GroupBy(_.d).Get(_.d, _.d_desc[0], _.genegroup.Unique()).Sort(_.genegroup.Count(), descend=True);
KEGG_genes_P1    = KEGG_pathways_P1.genegroup.Flat().Unique()

KEGG_pathways_P2 = KEGG[_.genegroup.In(diff_reg_all_H2.id)].GroupBy(_.d).Get(_.d, _.d_desc[0], _.genegroup.Unique()).Sort(_.genegroup.Count(), descend=True);
KEGG_genes_P2    = KEGG_pathways_P2.genegroup.Flat().Unique()

PATHWAYS_AFFECTED = KEGG[_.genegroup.In(diff_reg_all.id)].GroupBy(_.d).Unique(_.genegroup).Get(_.d, _.d_desc[0], _.genegroup.Count(), _.genegroup[_.In(diff_reg_all_H1.id)] / 'h1', _.genegroup[_.In(diff_reg_all_H2.id)] / 'h2').Sort(_.genegroup, descend=True).Show()

Export(PATHWAYS_AFFECTED, '%s/kegg.affected_pathways.tsv' % output_dir)

###############################################################################
###############################################################################
###############################################################################
  # Statistical test for P1/P2 diff ex significance

def fdr_gen(p, alpha, type='bh'):

  nt     = len(p);  # Number of tests
  ps   = sorted(p);
  indx = [ i[0] for i in sorted(enumerate(p), key=lambda x:x[1]) ];

  if type == 'bhy':
    cm    = sum([ 1.0/float(i) for i in xrange(nt)] );
    klist = [ (float(i+1)/(float(nt) * cm)) for i in xrange(nt) ];
  else:
    klist = [ (float(i+1)/float(nt)) for i in xrange(nt) ];
  #fi

    # Adjust pvalues to qvalues
  padj = [ ps[i] / klist[i] for i in xrange(nt)];
    # Fix pvalues larger than 1
  q = [ qi if qi < 1.0 else 1.0 for qi in padj ];

    # Monotonicity
  qm = [];
  prev_v = q[0];
  for v in q:
    qm.append(max(prev_v, v));
    prev_v = qm[-1];
  #efor

    # get back to original sorting
  qrs = [0] * nt;
  for i in xrange(nt):
    qrs[indx[i]] = qm[i];
  #efor

  return qrs;
#edef

###############################################################################

from scipy.stats import binom
diff_reg_counts     = diff_reg_all.GroupBy(_.condition).Get(_.condition, _.foldchange[_.foldchange < 1.0/3.0].Count() / 'h1', _.foldchange[_.foldchange > 3].Count() / 'h2')
diff_reg_counts_all = diff_reg_counts | Stack | Rep([('mushroom', 176, 192), ('compost', 30, 52)])

h0_mu = 0.5

binom_tests = []

for (condition, h1_up, h2_up) in zip(*diff_reg_counts_all()):
  binom_tests.append((condition, h1_up, h2_up, 1-binom.cdf(max(h1_up,h2_up), h1_up+h2_up, h0_mu)))
#efor

binom_pvalues     = [ x[-1] for x in binom_tests ]
binom_fdr_pvalues = fdr_gen(binom_pvalues, 0.05)

corr_binom_tests = []
for (condition, h1_up, h2_up, pvalue), padj in zip(binom_tests, binom_fdr_pvalues):
  corr_binom_tests.append((condition, h1_up, h2_up, pvalue, padj))
#efor

corr_binom_tests_rep = Rep(corr_binom_tests) / ('condition', 'p1_up', 'p2_up', 'pvalue', 'qvalue')

Export(corr_binom_tests_rep.Sort(_.qvalue), '%s/diff.diff_reg_test.tsv' % output_dir)
#########################################################################################
#########################################################################################
#########################################################################################
  # Plot number of markers

np1, np2, S = (Mflat | Match(_.geneid, _.geneid) | D).GroupBy(_.genegroup).Sort(_.geneid).Unique(_.geneid)[_.nmarkers.Min() > 5].nmarkers[_.nmarkers.Count() == 2].Get(_.nmarkers[0] / 'np1', _.nmarkers[1] / 'np2', _.nmarkers[1].Each(lambda x: 'c').Cast(str)).GroupBy((_.np1, _.np2)).Get(_.np1[0], _.np2[0], _.nmarkers.Count())()

plt.clf(); plt.cla();

plt.figure(figsize=(8,8))
sns.set(style="white", color_codes=True)
g = plt.scatter(np1, np2, s=S)
plt.xlabel('Number of markers in P1 gene')
plt.ylabel('Number of markers in P2 gene');

plt.xscale("log", nonposx='clip')
plt.yscale("log", nonposy='clip')

minval = 0
maxval = max(max(np1), max(np2)) + 5

plt.xlim([minval, maxval])
plt.ylim([minval, maxval])

plt.savefig('%s/diff.nmarkers.png' % output_dir)
plt.savefig('%s/diff.nmarkers.svg' % output_dir)

#########################################################################################
#########################################################################################
#########################################################################################
  # Average ratios being different.

from scipy.stats import t as tdist;

ratios = DI_chr_ratio.mean.Each(lambda x: np.log(x)).Cast(float)()

ratios_mean = np.mean(ratios)
ratios_stdev = np.std(ratios);

t = (ratios_mean-0) / (ratios_stdev / np.sqrt(ratios.shape[0]))

1-tdist.cdf(t, ratios.shape[0]-1)

#########################################################################################
#########################################################################################
#########################################################################################
  # Is the number of markers correlated to expression?

from scipy.stats import linregress

nmarkers    = Dgenegroup.GroupBy(_.geneid).Get(_.geneid, _.nmarkers[0]);
max_expr    = (DT_genes.Get(_.id.Each(lambda x: x.split(',')[0]).Cast(str), _.basemeana) | Stack | DT_genes.Get(_.id.Each(lambda x: x.split(',')[1]).Cast(str), _.basemeanb)).GroupBy(_.id).Get(_.id, _.result.Max() / 'expr')
marker_expr = max_expr | Match(_.id, _.geneid, merge_same='equi') | nmarkers

plt.cla(); plt.clf();

X_nmarkers = marker_expr.nmarkers();
Y_expr     = marker_expr.expr()

linregress(X_nmarkers, Y_expr)

sns.regplot(x=X_nmarkers, y=Y_expr);

plt.xlabel('Number of markers');
plt.ylabel('Mean P1/P2 expression');

plt.savefig('%s/dist.nmarkers_expr.png' % output_dir)
plt.savefig('%s/dist.nmarkers_expr.svg' % output_dir)

#########################################################################################
#########################################################################################
#########################################################################################


