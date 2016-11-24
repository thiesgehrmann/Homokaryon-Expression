
import genome_drawing_split_meth as gdsm;
import python_utils as utils

import pandas as pd
import matplotlib.pylab as plt;
import seaborn as sns;

import numpy as np

from bx.intervals.intersection import IntervalTree

genome_chromosomes = Rep([('scaffold_1', 1, 3549961),
                          ('scaffold_10', 1, 1708347),
                          ('scaffold_11', 1, 1688363),
                          ('scaffold_12', 1, 1482390),
                          ('scaffold_13', 1, 1333233),
                          ('scaffold_2', 1, 3489538),
                          ('scaffold_3', 1, 3131717),
                          ('scaffold_4', 1, 3112497),
                          ('scaffold_5', 1, 2550519),
                          ('scaffold_6', 1, 2329483),
                          ('scaffold_7', 1, 2323014),
                          ('scaffold_8', 1, 1953080),
                          ('scaffold_9', 1, 1762519)]) / ('f0', 'start', 'end');

output_dir = './A15'
#METH = Read('/home/nfs/thiesgehrmann/groups/w/phd/data/methylation/A15/A15_CpG_methylated_H97v3.1.gff', format='tsv').Cast(str, str, str, int, int, int, str, str, str)[_.f0 != ''];
METH = Read('/home/nfs/thiesgehrmann/groups/w/phd/tasks/methylation/A15/methylated_bases.tsv', format='tsv').Cast(str, int, int, str, int, int, int) / ('seqname', 'start', 'end', 'strand', 'coverage', 'numc', 'numt');
METH = METH.To(_.seqname, Do=_.Each(lambda x: '_'.join(x.split('_')[:-1])))
METH = METH[_.numc > 0]

plt.cla(); plt.clf();
ax = sns.distplot( METH[_.numc != _.coverage].Get((_.numc.Cast(float) / _.coverage)).Cast(float)(), kde=False, bins=100);
#ax.set_yscale("log", nonposy='clip')
plt.xlabel("Fraction of methylated Cytosines vs. nonmethylated cytosines");
plt.ylabel("Frequency");
plt.savefig('/home/nfs/thiesgehrmann/groups/w/phd/tasks/methylation/A15/methylation_VAF.png')
plt.savefig('/home/nfs/thiesgehrmann/groups/w/phd/tasks/methylation/A15/methylation_VAF.svg')

DIFFMETH = METH.Get((_.numt.Cast(float) / _.coverage) / 'frac' , *METH.Names)
DIFFMETH = DIFFMETH[_.frac >= 0.4][_.frac < 0.6].Copy()
METH = METH.Copy()

COMPOST_genes = Read('A15_compost_nodup/gene_ratio_stats.tsv') / ('genegroup', 'cond', 'r1', 'r2', 'mean', 'sd', 'var', 'confidence', 'seqname', 'start', 'end')
COMPOST_genes = COMPOST_genes.Cast(str, str, float, float, float, float, float, float, str, int, int);
JORDI_genes   = Read('A15/gene_ratio_stats.tsv') / ('genegroup', 'cond', 'r1', 'r2', 'mean', 'sd', 'var', 'confidence', 'seqname', 'start', 'end');
JORDI_genes   = JORDI_genes.Cast(str, str, float, float, float, float, float, float, str, int, int);

COMPOST_DIFFREG = Read('A15_compost_nodup/diff_regulated.tsv').Detect();
JORDI_DIFFREG   = Read('A15/diff_regulated.tsv').Detect()

DIFFREG = COMPOST_DIFFREG | Stack | JORDI_DIFFREG

genes = COMPOST_genes | Stack | JORDI_genes

ideogram = genome_chromosomes.Get(_.f0, _.start, _.end, _.f0 / 'name', _.f0.Each(lambda x: (1,1,1)).Detect() / 'color');
ideogene = genes.GroupBy(_.genegroup).Get(_.genegroup, _.seqname[0], _.start[0], _.end[0], _.Get(_.mean.Min(), _.mean.Max()).Each(lambda x,y: x if np.abs(np.log(x)) > np.abs(np.log(y)) else y ).Cast(float) / 'mean');
ideogene = ideogene.Get(_.seqname, _.start, _.end, _.genegroup, _.mean,
                                _.mean.Each(lambda x: utils.rgb2frac(utils.convert_to_rgb(0.0, 1.0, utils.norm(0.5, 1.0, 2.0, x)))) / 'colors')
#ideometh = METH.Get(_.f0, _.f3) / ('seqname', 'start');
ideometh = METH.Get(_.seqname, _.start).Copy()

gdsm.draw_genome(utils.natural_sort(genome_chromosomes.f0.Unique()()),
               utils.rep_to_df(ideogram),
               utils.rep_to_df(ideogene),
               utils.rep_to_df(ideometh),
               '%s/meth.genes' % (output_dir))

###############################################################################
###############################################################################
###############################################################################

  # See which genes are near methylated bases

from scipy.stats import chi2_contingency

def index_genes(G, window=0):

  G = G.GroupBy(_.seqname).Sort(_.start);
  G = G.Get(_.seqname, _.name, _.start, _.end).Flat();

  chrs = {};
  for (seqname, name, start, end) in zip(*G()):
    if seqname not in chrs:
      chrs[seqname] = IntervalTree();
    #fi

    chrs[seqname].add(start-window, end+window, (name, start, end));
    print seqname, start, end, name;
  #efor

  return chrs;
#edef

###############################################################################

index = index_genes(genes.Unique(_.genegroup).Get(_.seqname, _.start, _.end, _.genegroup / 'name'), window=1000)

def get_meth_genes(METH, index):

  meth_genes = [];

  for meth_seq, meth_start in zip(*METH.Get(_.seqname, _.start)()):
    ov_genes = index[meth_seq].find(meth_start, meth_start);
    for g in ov_genes:
      meth_genes.append((meth_seq, meth_start, g[0], g[1], g[2]));
    #efor
  #efor

  MGENES = Rep(meth_genes) / ('seqname', 'base', 'id', 'start', 'end')
  MGENES = MGENES.GroupBy(_.id).Get(_.id, _.seqname[0], _.start[0], _.end[0], _.seqname.Count() / 'count').Sort(_.count, descend=True)
  return MGENES
#edef

MGENES     = get_meth_genes(METH, index);
DIFFMGENES = get_meth_genes(DIFFMETH, index)

sigdiffreg_all      = set(DIFFREG[_.padj < 0.05].id.Unique()())
sigdiffreg_compost  = set(COMPOST_DIFFREG[_.padj < 0.05].id.Unique()())
sigdiffreg_jordi    = set(JORDI_DIFFREG[_.padj < 0.05].id.Unique()())
sigdiffreg_overlap  = sigdiffreg_compost & sigdiffreg_jordi
sifdiffreg_com_uniq = sigdiffreg_compost - sigdiffreg_overlap
sigdiffreg_jor_uniq = sigdiffreg_jordi - sigdiffreg_overlap

def enrich_methylation(genes, sigdiffreg, MGENES):

  totalgenes   = set(genes.genegroup.Unique()())
  totaldiff    = set(sigdiffreg);
  totaldiff_no = totalgenes - totaldiff
  totalmeth    = set(MGENES.id());
  totalmeth_no = totalgenes - totalmeth

  diff_meth     = totaldiff & totalmeth
  diff_nometh   = totaldiff & totalmeth_no
  nodiff_meth   = totaldiff_no & totalmeth
  nodiff_nometh = totaldiff_no & totalmeth_no

  return chi2_contingency([ [len(diff_meth), len(diff_nometh)], [len(nodiff_meth), len(nodiff_nometh)] ]), [len(diff_meth), len(diff_nometh)], [len(nodiff_meth), len(nodiff_nometh)]
#edef

enrichments_mgenes = [ enrich_methylation(genes, x, MGENES) for x in [ sigdiffreg_all, sigdiffreg_compost, sigdiffreg_jordi, sigdiffreg_overlap, sifdiffreg_com_uniq, sigdiffreg_jor_uniq]]
[ [x[0][1]] + x[1] + x[2] for x in enrichments]

enrichments_diffmgenes = [ enrich_methylation(genes, x, DIFFMGENES) for x in [ sigdiffreg_all, sigdiffreg_compost, sigdiffreg_jordi, sigdiffreg_overlap, sifdiffreg_com_uniq, sigdiffreg_jor_uniq]]
[ [x[0][1]] + x[1] + x[2] for x in enrichments_diffmgenes]

