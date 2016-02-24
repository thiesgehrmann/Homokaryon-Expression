
import numpy as np
from scipy.stats import t

import python_utils as utils;
import python_svg_utils as svg_utils;
import pandas as pd
import matplotlib.pylab as plt;
import seaborn as sns;
###############################################################################

# The structure of the data

#sys.argv = ['',
#            'agabi_tissue_data.tsv',
#            "/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.assembly.fasta",
#            "/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.genes.gff3",
#            'A15',
#            'homokaryon_activity.svg']

sys.argv = [ '',
             'agabi_compost_data.tsv',
             "/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.assembly.fasta",
             "/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.genes.gff3",
             'A15_compost'];

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

###############################################################################
# Based on the mapping, create a view where we see both 

M = Read('%s/mapping.tsv' % (output_dir))

Mgroup = M.Get(_.Get(tuple(M.Names)).To(0, Do=_.Each(lambda x: ','.join(x)).Cast(str) / 'group'), *M.Names);
Mflat = Mgroup.Get(_.group, M.Names[0])
for orgname in M.Names[1:]:
  Mflat = Mflat | Stack | Mgroup.Get(_.group, orgname);
#efor
Mflat = Mflat / ('genegroup', 'geneid');

# Groups with genes 
Dgenegroup = Mflat | Match(_.geneid, _.geneid) | D
Dgenegroup = Dgenegroup.Sort(_.geneid, _.replicate_id);

fieldnames = Dgenegroup.GroupBy(_.genegroup).Sort(_.geneid, _.replicate_id)[0].Get(_.geneid, _.sample_name).To(_.geneid, Do=_.Each(lambda x: x.split('|')[0]).Cast(str)).Get(_.geneid + '|' + _.sample_name / 'fieldnames')().tolist()
fieldnames = tuple([ x.lower() for x in ['genegroup'] + fieldnames ])
deseq_input = Rep([ (x,) + tuple(y.tolist()) for (x,y) in zip(*Dgenegroup.Get(_.genegroup, _.avg_read_depth.Cast(int)).GroupBy(_.genegroup)())]) / fieldnames
Export(deseq_input, '%s/DESEQ_input.tsv' % (output_dir));

###############################################################################
# HERE WE NEED TO RUN DESEQ/EdgeR/whatever for normalization!!!
cmd = "Rscript deseq.R '%s/DESEQ_input.tsv' %s" % (output_dir, output_dir);
utils.run_cmd(cmd)

###############################################################################

DT = Read('%s/deseq_tests.tsv' % output_dir).Detect()

# We have two homokaryons, Homokaryon1 and Homokaryon2, H1 and H2

DT_H1 = (DT.To(_.id, Do=_.Each(lambda x: x.split(',')[0]).Cast(str)) | Match(_.id, jointype='left', merge_same=True) | gffs[0].Get(_.seqname, _.id))
DT_H2 = (DT.To(_.id, Do=_.Each(lambda x: x.split(',')[1]).Cast(str)) | Match(_.id, jointype='left', merge_same=True) | gffs[1].Get(_.seqname, _.id))

all_counted_ids = (DT_H1.id | Stack | DT_H2.id).Copy();

upregulated_in_H1 = DT_H1[_.padj < .05][_.foldchange < 1.0/3.0]
upregulated_in_H2 = DT_H2[_.padj < .05][_.foldchange > 3]

def format_upregulated_conds(data):
  conds = data.GroupBy(_.condition).Get(_.condition, _.id, _.seqname)
  any   = data.To(_.condition, Do=_.Each(lambda x: 'any').Cast(str)).Get(_.condition, _.id, _.seqname).Unique().GroupBy(_.condition)
  all   = data.GroupBy(_.id).Get(_.id, _.seqname.Count() / 'n', _.seqname[0])[_.n == len(label_names)].GroupBy(_.n).Get(_.n.Each(lambda x: 'all').Cast(str), _.id, _.seqname)
  return ((conds | Stack | any ) | Stack | all) / ('cond', 'id', 'seqname');
#edef

from scipy import stats as ssp;

def enrich(data):
  data = data / ('name', 'a', 'b', 'c', 'd')
  tested = data.Get(_.name, _.a, _.b, _.c, _.d, _.Get(_.a, _.b, _.c, _.d).Each(lambda a,b,c,d: ssp.fisher_exact([ [a, max(b, 1)], [max(c,1), d] ])) / 'test')
  return tested.Get(_.name, _.a, _.b, _.c, _.d, _.test.Each(lambda x: x[0]).Cast(float) / 'r', _.test.Each(lambda x: x[1]).Cast(float) / 'p').Copy()
#edef

def enriched_chromosome(data):
  data = data / ('seqname', 'ndiff', 'ngenes')
  total_diff  = data.ndiff.Sum()()
  total_genes = data.ngenes.Sum()()
  enrich_data = data.Get(_.seqname, _.ndiff / 'a', (_.ngenes - _.ndiff) / 'b', _.ndiff.Each(lambda x: total_diff - x).Cast(int) / 'c', _.Get(_.ngenes, _.ndiff).Each(lambda x, y: (total_genes - (x - y))) / 'd')
  return enrich(enrich_data)
#edef

upregulated_conds_H1 = format_upregulated_conds(upregulated_in_H1)
upregulated_conds_H2 = format_upregulated_conds(upregulated_in_H2)

all_enrichtests = []
for org, upregulated_conds, genome in zip(['H1', 'H2'], [ upregulated_conds_H1, upregulated_conds_H2], gffs):
  for cond in upregulated_conds.cond():
    print cond
    data_left  = upregulated_conds.Flat()[_.cond == cond].GroupBy((_.cond, _.seqname)).Get(_.cond[0], _.id.Count() / 'ndiff', _.seqname[0])
    data_right = genome.GroupBy(_.seqname)[_.feature == 'gene'][_.id.In(all_counted_ids)].Get(_.seqname, _.id.Count() / 'ngenes')
    data = data_left | Match(_.seqname, jointype='full', merge_same=True) | data_right
    data = data.Get(_.seqname, _.ndiff, _.ngenes).ReplaceMissing()
    all_enrichtests.append(enriched_chromosome(data).To(_.name, Do=_.Each(lambda x: org + '|' + cond + '|' + x).Cast(str)).Copy())
  #efor
#efor

enrichtests = all_enrichtests[0]
for et in all_enrichtests[1:]:
  enrichtests = enrichtests | Stack | et;
#efor

enrichtests[_.p < 0.05 / enrichtests.Shape()()].Show()

scaffold12b_genes = upregulated_conds_H1[_.seqname == 'scaffold_12b'].id.Flat().Unique().Sort().Show()
mapping0 = Read('%s/mapping_0.fasta' % output_dir)

###############################################################################





mapping0[_.f0.In(scaffold12b_genes.id)].Sort(_.f0).To(_.seq, Do=_.Each(lambda x: utils.translate(x)).Cast('protein'))

###############################################################################

DI_chr = deseq_input | Match(_.genegroup.Each(lambda x: x.split(',')[0]).Cast(str), _.id, merge_same='equi', jointype='left') | gffs[0].Get(_.id, _.seqname)
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
  cond_stats = cond_stats.Get(_.seqname, _.r1, _.r2, ((_.r1 + _.r2) / 2) / 'mean' )
  cond_stats = cond_stats.Get(_.seqname, _.r1, _.r2, _.mean, ((_.r1 - _.mean) * (_.r1 - _.mean) + (_.r2 - _.mean) * (_.r2 - _.mean) / 2) / 'var')
  cond_stats = cond_stats.Get(_.seqname, _.r1, _.r2, _.mean, _.var.Each(lambda x: np.sqrt(x)).Cast(float) / 'sd', _.var );
  cond_stats = cond_stats.Get(_.seqname, _.r1, _.r2, _.mean, _.sd, _.var, (_.sd / np.sqrt(2)) * max(t.interval(0.90, 1)) / 'confidence' ).Copy()
  stats.append(cond_stats.Get(_.seqname.Each(lambda x: cond).Cast(str) / 'cond', *cond_stats.Names).Copy());
#efor

allstats = stats[0];
for s in stats[1:]:
  allstats = allstats | Stack | s;
#efor

  # Add the data from ALL the chromosomes
condlibratios = Read('%s/deseq_condlibratios.tsv' % output_dir).Detect() / ('cond', 'r1', 'r2', 'mean', 'sd', 'var');
condlibratios = condlibratios.Get(_.cond, _.r1, _.r2, _.mean, _.sd, _.var, (_.sd / np.sqrt(2)) * max(t.interval(0.90, 1)) / 'confidence')

allstats = condlibratios.Get(_.cond, _.cond.Each(lambda x: 'all').Cast(str) / 'seqname', _.r1, _.r2, _.mean, _.sd, _.var, _.confidence) | Stack | allstats

minval, maxval = allstats.Get(_.mean.Min(), _.mean.Max())()
colors = [(0, 0, 255), (255, 255, 255), (255, 0, 0)]  # [BLUE, WHITE, RED]

defaults = { '!!H1!!': 'P1', '!!H2!!': 'P2', '!!ratio_min!!': '%1.1f' % minval, '!!ratio_max!!': '%1.1f' % maxval }

def draw_heatmap(allstats, outfigure, order=None):
  conds, seqname, mean = allstats.Get(_.cond, _.seqname, _.mean)();
  df = pd.DataFrame({ 'conds': conds, 'seqname': seqname, 'mean':mean})
  df = df.pivot('conds', 'seqname', 'mean')
  plt.cla(); plt.clf();
  ax = sns.heatmap(df, center=1.0)
  plt.savefig('%s/condlibratios.svg' % output_dir)
#edef

draw_heatmap(allstats, '%s/condlibratios.svg' % output_dir)

# Scale so that 1.0 is in the middle...
#We need to scale in a piece wise continuous way...
def norm(min, mid, max, x):

  nmin = 0;
  nmid = 0.5
  nmax = 1.0

  if x < mid:
    return  (x-min) * (nmid / (mid - min))
  else:
    return (x-mid) * ((nmax - nmid) / (max - mid)) + nmid
  #fi
#edef

  


allstats_mean_norm = allstats.Get(_.mean.Each(lambda x: norm(minval, 1.0, maxval, x)).Cast(float) / 'norm', *allstats.Names)

if tissue_figure is not None:

  reload(svg_utils)
  reload(utils);
  for chrid in allstats.seqname.Unique()():
    replacevals = defaults
    replacevals.update(dict(zip(*allstats_mean_norm[_.seqname == chrid].Get(_.cond.Each(lambda x: '!!' + x + '!!').Cast(str), 
                                                                            _.norm.Each(lambda x: utils.rgb2hex(utils.convert_to_rgb(0, 1, x, colors))))())))
    svg = svg_utils.read_svg(tissue_figure);
    svg = svg_utils.replacem(svg, replacevals);
    svg_utils.write_svg(svg, '%s/%s.homokaryon_activity.svg' % (output_dir, chrid))
  #efor
  
  
  for chrid in allstats.seqname.Unique()():
    replacevals = defaults
    replacevals.update(dict(zip(*allstats_mean_norm[_.seqname == chrid].Get(_.cond.Each(lambda x: '!!' + x + '!!').Cast(str),
                                                                            _.Get(_.norm, _.mean, _.confidence).Each(lambda w,x,y: utils.rgb2hex(utils.convert_to_rgb(0, 1, w, colors)) if np.abs(1-x) > y else '000000' ))())))
    print chrid
    print allstats_mean_norm[_.seqname == chrid].Get(_.cond, _.Get(_.norm, _.mean, _.confidence).Each(lambda w,x,y: (w,x,y, utils.rgb2hex(utils.convert_to_rgb(0, 1, w, colors)) if np.abs(1-x) > y else '000000' )))()
    svg = svg_utils.read_svg(tissue_figure);
    svg = svg_utils.replacem(svg, replacevals);
    svg_utils.write_svg(svg, '%s/%s.signif.homokaryon_activity.svg' % (output_dir, chrid))
  #efor
#fi

###############################################################################
