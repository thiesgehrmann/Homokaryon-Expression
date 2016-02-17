
import numpy as np

###############################################################################

# The structure of the data

sample_names  = np.array([ "Casing_C_26" ,"Initials_I272" ,"Pileal_Stipeal_A2611" ,"Pileal_Stipeal_A2612" ,"Different_D2641" ,"Different_D2642" ,"Different_D2643" ,"Different_D2644" ,"Different_D2645" ,"Different_D2646" ,"YFB_F2631" ,"YFB_F2632" ,"YFB_F2633" ,"YFB_F2634" ,"YFB_F2635" ,"YFB_F2636" ,"YFB_F2637" ,"Casing_C_27" ,"Initials_I2611" ,"Pileal_Stipeal_A2711" ,"Pileal_Stipeal_A2712" ,"Different_D2721" ,"Different_D2722" ,"Different_D2723" ,"Different_D2724" ,"Different_D2725" ,"Different_D2726" ,"YFB_F2731" ,"YFB_F2732" ,"YFB_F2733" ,"YFB_F2734" ,"YFB_F2735" ,"YFB_F2736" ,"YFB_F2737"]);
  # Replicate groups
sample_labels = np.array([ 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 ]);

sample_names = np.concatenate([ sample_names[sample_labels == i] for i in xrange(max(sample_labels+1))])
sample_labels = np.concatenate([ sample_labels[sample_labels == i] for i in xrange(max(sample_labels+1))])

  # Replicate group names
label_names   = np.array([ "Veg","Initials","PS_stipe_center","PS_stipe_shell","DIF_stipe_center","DIF_stipe_shell","DIF_stipe_skin","DIF_cap_skin","DIF_cap_tissue","DIF_gill_tissue","YFB_stipe_center","YFB_stipe_shell","YFB_stipe_skin","YFB_cap_skin","YFB_cap_tissue","YFB_gill_tissue","YFB_veil"]);

###############################################################################

# Gather all the read counts, calculate average per gene, and put into a big table

allD = [];
allDs = []
for i, label_name in enumerate(label_names):
  samples    = sample_names[sample_labels == i]
  sample_ids = np.where(sample_labels == i)[0]
  print samples;
  Ds = [ Read("FIXTHIS_%s.tsv.ara" % (x), verbose=False).Detect() for x in samples];
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

M = Read('mapping.tsv')

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
Export(deseq_input, 'DESEQ_input.tsv');

# HERE WE NEED TO RUN DESEQ/EdgeR/whatever for normalization!!!

Dgenegroup.GroupBy((_.geneid, _.sample_name))

###############################################################################

DT = Read('deseq_tests.tsv').Detect()

genomes = [ Read(x, format='gff') for x in '/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h39_1/AgabiH39_1.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h97_1/AgabiH97_1.gff3'.split(',')]


DT_H39 = (DT.To(_.id, Do=_.Each(lambda x: x.split(',')[0]).Cast(str)) | Match(_.id, jointype='left', merge_same=True) | genomes[0].Get(_.seqname, _.id))
DT_H97 = (DT.To(_.id, Do=_.Each(lambda x: x.split(',')[1]).Cast(str)) | Match(_.id, jointype='left', merge_same=True) | genomes[1].Get(_.seqname, _.id))

upregulated_in_h39 = DT_H39[_.padj < .05][_.foldchange < 1.0/3.0]
upregulated_in_h97 = DT_H97[_.padj < .05][_.foldchange > 3]

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

upregulated_conds_h39 = format_upregulated_conds(upregulated_in_h39)
upregulated_conds_h97 = format_upregulated_conds(upregulated_in_h97)

all_enrichtests = []
for org, upregulated_conds, genome in zip(['h39', 'h97'], [ upregulated_conds_h39, upregulated_conds_h97], genomes):
  for cond in upregulated_conds.cond():
    print cond
    data = upregulated_conds.Flat()[_.cond == cond].GroupBy((_.cond, _.seqname)).Get(_.cond[0], _.id.Count() / 'ndiff', _.seqname[0]) | Match(_.seqname, jointype='full', merge_same=True) | genome.GroupBy(_.seqname)[_.feature == 'gene'].Get(_.seqname, _.id.Count() / 'ngenes')
    data = data.Get(_.seqname, _.ndiff, _.ngenes).ReplaceMissing()
    all_enrichtests.append(enriched_chromosome(data).To(_.name, Do=_.Each(lambda x: org + '|' + cond + '|' + x).Cast(str)).Copy())
  #efor
#efor

enrichtests = all_enrichtests[0]
for et in all_enrichtests[1:]:
  enrichtests = enrichtests | Stack | et;
#efor

enrichtests[_.p < 0.05 / enrichtests.Shape()()].Show()

scaffold12b_genes = upregulated_conds_h39[_.seqname == 'scaffold_12b'].id.Flat().Unique().Sort().Show()
genes = Read('mapping_0.fasta')

bases = ['t', 'c', 'a', 'g'];
codons = [a+b+c for a in bases for b in bases for c in bases];
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
codon_table = dict(zip(codons, amino_acids));

def translate(seq):
  aa = '';
  for i in xrange(0, len(seq), 3):
    cod = seq[i:i+3].lower();
    aa += codon_table[cod] if (cod in codon_table) else '*';
  #efor
  return aa;
#edef

genes[_.f0.In(scaffold12b_genes.id)].Sort(_.f0).To(_.seq, Do=_.Each(lambda x: translate(x)).Cast('protein'))

###############################################################################

condlibratios = Read('deseq_condlibratios.tsv').Detect()

from scipy.stats import t

condlibratios = condlibratios.Get(_.cond, _.r1, _.r2, _.avg, _.sd, _.var, (_.sd / np.sqrt(2)) * max(t.interval(0.95, 1)) / 'confidence')

# from http://stackoverflow.com/questions/20792445/calculate-rgb-value-for-a-range-of-values-to-create-heat-map
def convert_to_rgb(minval, maxval, val, colors):
    max_index = len(colors)-1
    v = float(val-minval) / float(maxval-minval) * max_index
    i1, i2 = int(v), min(int(v)+1, max_index)
    (r1, g1, b1), (r2, g2, b2) = colors[i1], colors[i2]
    f = v - i1
    return int(r1 + f*(r2-r1)), int(g1 + f*(g2-g1)), int(b1 + f*(b2-b1))

minval, maxval = condlibratios.avg.Min()(), condlibratios.avg.Max()()
colors = [(0, 0, 255), (255, 255, 255), (255, 0, 0)]  # [BLUE, GREEN, RED]
condlibratios.Get(_.cond, _.avg.Each(lambda x: convert_to_rgb(minval, maxval, x, colors)))




