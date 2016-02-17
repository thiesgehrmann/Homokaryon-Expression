#!/bin/env python
####################################################################
# Create a mapping between
####################################################################


from ibidas import *;
from external_packages.proteny import proteny as ps;
from external_packages.proteny import cluster_null;
from external_packages.proteny import data as pdata;
import sys;
import os;

sys.path.append('./external_packages/proteny')
import cluster_null as cluster_null;

####################################################################
  # Usage function

def usage(arg0):
  print "Mapping Usage:";
  print " %s <org_genome_files> <org_gff_files> <output_dir>" % arg0;
  print "";
  print "org_genome_files: (files) Comma-separated list of FASTA files with DNA genome sequences"
  print "org_gff_files:    (files) Comma-separated list of GFF3 files with gene annotations for the genome sequences given in <org_genome_files>"
  print "";
#edef

####################################################################
  # Prepare data for proteny

def prepare_proteny(name, gff, seq):
  S = Read(seq, format='fasta');
  G = Read(gff, format='gff');

  GE = G[_.feature == 'exon'];
  GE = GE.To(_.id, Do=_.Each(lambda x: x.split('|')[-1][4:]).Cast(int));
  GE = GE.To(_.parent, Do=_.Each(lambda x: x.split('|')[-1]).Cast(str));
  GE = GE.Get(_.seqname, _.start, _.end, _.strand, _.parent, _.parent, _.id).Copy();

  return (name, GE / pdata.__gene_slice_names__, S / pdata.__genome_slice_names__);
#edef

####################################################################
  # Prepare for DNA blast

from string import maketrans;
def revcomp(dna_seq):
  return dna_seq[::-1].translate(maketrans("ATGC","TACG"))
#edef

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

def prot_seqs_from_gff(G, S):
  GS = G[_.feature == 'CDS'] | Match(_.seqname, 0) | S
  GS = GS.To(_.seq, Do=_.Get(GS.start, GS.end, GS.seq).Each(lambda w,x,y: y[w-1:x]) / 'seq');
  GS = GS.GroupBy(_.parent).Get(_.seqname[0], _.source[0], _.id[0], _.parent, _.name[0], _.feature[0], _.start.Min(), _.end.Max(), _.score[0], _.strand[0], _.frame[0], _.attr, _.seq.Sum());
  GS = GS.To(_.seq, Do=_.Get(GS.strand, GS.seq).Each(lambda x,y: translate(revcomp(y) if x == '-' else y)) / 'seq');
  GS = GS.To(_.seq, Do=_.Cast('protein'));
  return GS.Get(_.parent, _.seq);
#edef


def dna_seqs_from_gff(G, S):
  GS = G[_.feature == 'CDS'] | Match(_.seqname, 0) | S 
  GS = GS.To(_.seq, Do=_.Get(GS.start, GS.end, GS.seq).Each(lambda w,x,y: y[w-1:x]) / 'seq');
  GS = GS.GroupBy(_.parent).Get(_.seqname[0], _.source[0], _.id[0], _.parent, _.name[0], _.feature[0], _.start.Min(), _.end.Max(), _.score[0], _.strand[0], _.frame[0], _.attr, _.seq.Sum());
  GS = GS.To(_.seq, Do=_.Get(GS.strand, GS.seq).Each(lambda x,y: revcomp(y) if x == '-' else y) / 'seq');
  GS = GS.To(_.seq, Do=_.Cast('DNA'));
  return GS.Get(_.parent, _.seq);
#edef

def blast_to_dictindex(BR):

  hits = zip(*BR());

  Hdict = {}

  for h in hits:
    if h[0] in Hdict:
      Hdict[h[0]] = Hdict[h[0]] + [ h ];
    else:
      Hdict[h[0]] = [h];
    #fi
  #efor

  Hdict = { k : sorted(Hdict[k], key=lambda x: x[12]) for k in Hdict };

  return Hdict;

#edef

####################################################################

def reciprocal_blast(BR_ab, BR_ba):
  d_ab = blast_to_dictindex(BR_ab)
  d_ba = blast_to_dictindex(BR_ba)

  rBR = []

  for k in d_ab:
    k_tophit = d_ab[k][0][14]
    if (k_tophit in d_ba) and (d_ba[k_tophit][0][14] == k):
      rBR.append(d_ab[k][0]);
    #fi
  #efor
  return Rep(rBR[0:]) / ( ('parentl', 'seql') + tuple(BR_ab.Names[2:-2]) + ('parentr', 'seqr'))
#edef

####################################################################
####################################################################
####################################################################
####################################################################

#os.sys.argv = [ 'mapping.py', 
#                '/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h39_1/AgabiH39_1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h97_1/agabiH97_1.assembly.fasta',
#                '/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h39_1/AgabiH39_1.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/agabi_h97_1/AgabiH97_1.gff3',
#                'U1'];

sys.argv = [ sys.argv[0],
             '/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.assembly.fasta,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.assembly.fasta',
             '/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p1.genes.gff3,/home/nfs/thiesgehrmann/groups/w/phd/data/A15_homokaryons/AgabiA15p2.genes.gff3',
             '/home/nfs/thiesgehrmann/groups/w/phd/tasks/karyon_specific_expression/Homokaryon-Expression/A15']

# Check if input is ok
if len(sys.argv) != 4:
  usage(sys.argv[0]);
  sys.exit(1);
#fi

  # Read input data
genome_files   = sys.argv[1].split(',');
gff_files      = sys.argv[2].split(',');
output_dir     = sys.argv[3];
blast_dir      = '/home/nfs/thiesgehrmann/groups/w/phd/tasks/blast_runs';

print genome_files
print gff_files
print output_dir

####################################################################
#Blast the DNA sequences

org_dnas = [];

for org_genome, org_gff in zip(genome_files, gff_files):
  print org_genome, org_gff;
  org_genome = Read(org_genome, format='fasta');
  org_gff    = Read(org_gff, format='gff');
  org_dnas.append(dna_seqs_from_gff(org_gff, org_genome).Copy())
#efor

org_blasts = [];
for i in xrange(len(org_dnas)):
  for j in xrange(i+1, len(org_dnas)):
    BR_ab = (org_dnas[i] | Blast(reciprocal=False, normalize=False, overwrite=False, folder=blast_dir) | org_dnas[j]).Copy()
    BR_ba = (org_dnas[j] | Blast(reciprocal=False, normalize=False, overwrite=False, folder=blast_dir) | org_dnas[i]).Copy()
    rBR = reciprocal_blast(BR_ab, BR_ba);
    org_blasts.append(rBR)
  #efor
#efor

###
# Need to do some magic here if we want more than two species...
# However, if we only want two, then we do this:
###
unique_mapped_diff = org_blasts[0][_.evalue < 1e-100][_.pident != 100];
MBR = [ unique_mapped_diff.Get(_.parentl / 'parent', _.seql / 'seq'), unique_mapped_diff.Get(_.parentr / 'parent', _.seqr / 'seq') ]

Export(Rep(zip(*[ x.parent() for x in MBR ])), '%s/mapping.tsv' % (output_dir), names=False);

for i, M in enumerate(MBR):
  Export(M, '%s/mapping_%d.fasta' % (output_dir, i));
#efor

####################################################################
#Perform a proteny synteny analysis

#pvalue  = 0.05;
#cthresh = 3;

#PR              = ps.proteny();
#org_proteny_ids = [];
#for i, (genome_file, gff_file) in enumerate(zip(genome_files, gff_files)):
#  org_proteny_id = PR.add_org(*prepare_proteny('org_%d' % (i+1), gff_file, genome_file), isfile=False)
#  org_proteny_ids.append(org_proteny_id);
##efor
#k = PR.analyze(id_a=org_proteny_ids[0], id_b=org_proteny_ids[1], cut='deeper_greater', nd=cluster_null.cluster_null_score_strict_smart, alpha=pvalue, ngenes_threshold=2, conservation_ratio=cthresh);

