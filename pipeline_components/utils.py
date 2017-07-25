###############################################################################

import errno
import os
import csv
import gzip

# https://stackoverflow.com/a/600612
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

###############################################################################

from collections import namedtuple


class BlastHitType(object):
  BlastHit = None
  def __init__(self):
    self.BlastHit = namedtuple("BlastHit", "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore")
  def setFields(self, blastFields):
    self.BlastHit = namedtuple("BlastHit", blastFields)
  #edef
  def getFields(self):
    return self.BlastHit._fields

  def blastFieldTypeConversion(self, fieldName, value):
    typeMap = { "qseqid" : lambda x: str(x),
                "sseqid" : lambda x: str(x),
                "pident" : lambda x: float(x),
                "length" : lambda x: int(x),
                "mismatch" : lambda x: int(x),
                "gapopen": lambda x: int(x),
                "qstart": lambda x: int(x),
                "qend": lambda x: int(x),
                "sstart": lambda x: int(x),
                "send": lambda x: int(x),
                "evalue": lambda x: float(x),
                "bitscore": lambda x: float(x),
                "slen": lambda x: int(x),
                "qlen": lambda x: int(x) }
    return typeMap[fieldName](value)
  #edef
#eclass

blastHitType = BlastHitType()

def parseBlastHit(row):

  typedRow = [ blastHitType.blastFieldTypeConversion(field, value) for (field, value) in zip(blastHitType.getFields(), row) ]
  return blastHitType.BlastHit(*typedRow)
#edef

###############################################################################

def blastHit2Row(h):
  return [str(x) for x in [h[i] for i in range(len(h._fields))]]
#edef

###############################################################################

def readBlastFile(filename, fields="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"):
  blastHitType.setFields(fields)
  hits    = []
  # Read the BLAST hits
  with open(filename, "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t", quotechar="\"");
    for row in reader:
      hits.append(parseBlastHit(row))
    #efor
  #ewith
  return hits
#edef

###############################################################################

def writeBlastFile(H, filename):
  with open(filename, "w") as ofd_hits:
    for h in H:
      ofd_hits.write('\t'.join(blastHit2Row(h)) + '\n')
    #efor
  #ewith
#edef

###############################################################################

def indexListBy(L, key=lambda x: x[0]):
  G = {}
  for item in L:
    k = key(item)
    if k not in G:
      G[k] = []
    #fi
    G[k].append(item)
  #efor
  return G
#edef

###############################################################################

def blast_hits_overlap(hit_i, hit_j):
  if (hit_i.qseqid != hit_j.qseqid) or (hit_i.sseqid != hit_j.sseqid):
    return False
  else:
    # Overlapping either in the QUERY OR IN THE SUBJECT!
    return (max(hit_i.sstart, hit_j.sstart) <= min(hit_i.send, hit_j.send)) or (max(hit_i.qstart, hit_j.qstart) <= min(hit_i.qend, hit_j.qend))
#edef

def regions_overlap(r1, r2):
  (s1, e1) = r1
  (s2, e2) = r2
  # <= do not allow same start as end end. (100, 200), (200, 300) NOT ALLOWED
  # <  allow same start as end. (100, 200), (200, 300) ALLOWED
  return max(0, (min(e1, e2) - max(s1, s2)))
#edef

###############################################################################

def blast_hits_overlap_group_helper(hit_i, hit_j, allowed_overlap=5):
  same_query   = hit_i.qseqid == hit_j.qseqid
  same_subject = hit_i.sseqid == hit_j.sseqid

  query_overlap   = regions_overlap((hit_i.qstart, hit_i.qend), (hit_j.qstart, hit_j.qend)) > allowed_overlap
  subject_overlap = regions_overlap((hit_i.sstart, hit_i.send), (hit_j.sstart, hit_j.send)) > allowed_overlap

  return (same_query and query_overlap) or (same_subject and subject_overlap)
#wedef

def blast_hits_overlap_group(hits):
  return any([ blast_hits_overlap_group_helper(hits[i],hits[j]) for i in range(len(hits)-1) for j in range(i+1,len(hits)) ])
#edef

###############################################################################

import gzip

def loadFasta(fastaFile):

  F = {'': 0}

  current_seq = ""
  with (gzip.open(fastaFile, "r") if fastaFile[-2:] == "gz" else open(fastaFile, "r")) as fd:
    for line in fd:
      line = line.strip()
      if len(line) == 0:
        continue
      #fi
      if line[0] == '>':
        current_seq = line[1:].split(' ')[0]
        F[current_seq] = ""
      else:
        F[current_seq] = F[current_seq] + line.strip()
      #fi
  #ewith
  F.pop("", None)
  return F
#edef

###############################################################################

def writeFasta(fasta, outFile, linelength=80):
  with open(outFile, "w") as ofd:
    for  (name, sequence) in fasta:
      ofd.write(">%s\n" % name)
      ofd.write("%s\n" % '\n'.join([sequence[i:i+linelength] for i in range(0, len(sequence), linelength)]))
    #efor
  #ewith
#edef

###############################################################################

import numpy as np
def region_covered(regions, size=None):
  if size is None:
    size = max([ r[1] for r in regions ])
  #fi
  R = np.zeros(size)
  for r in regions:
    R[r[0]-1:r[1]-1] = 1
  #efor
  return R
#efor

###############################################################################

def region_coverage(regions, size = None):
  if size is None:
    size = max([ r[1] for r in regions ])
  #fi
  R = np.zeros(size)
  for r in regions:
    R[r[0]-1:r[1]-1] += 1
  #efor
  return R
#edef

###############################################################################

gff3Entry = namedtuple("gff3Entry", "seqid, source, type, start, end, score, strand, phase, attr")

def parseGFF3entry(fields):

   (seqid, source, feature, start, end, score, strand, phase, attr) = fields
   attr = dict([tuple(x.strip().split('=')[0:2]) for x in attr.split(";") ])
   return gff3Entry(seqid, source, feature, int(start), int(end), score, strand, phase, attr)
#edef

###############################################################################

from intervaltree import Interval, IntervalTree
class GFF3(object):
  def __init__(self, entries):
    self.entries = entries
    self.seqids  = set([ e.seqid for e in self.entries])
    self.genes   = { e.attr["ID"] : e for e in self.entries if e.type.lower() == 'mrna' }
    self.interval = self.indexByInterval()
    #self.geneIndex = self.indexByTopLevel()
  #edef

  def indexByTopLevelID(self):
    topLevel = { e.attr["ID"]: [e] for e in self.entries if "Parent" not in e.attr }
    lowerLevel = { e.attr["ID"] : e.attr["Parent"] for e in self.entries if "Parent" in e.attr }
    for e in [ e for e in self.entries if "Parent" in e.attr ]:
      parent = e.attr["Parent"]
      while parent not in topLevel:
        if parent in lowerLevel:
          parent = lowerLevel[parent]
        else:
          print("%s NOT in lowerLevel!" % parent)
          break
        #fi
      #ewhile
      topLevel[parent].append(e)
    #efor
    return topLevel
  #edef

  def indexByInterval(self):
    t = { seqid: IntervalTree() for seqid in self.seqids }
    for e in [ e for e in self.entries if e.type.lower() == "mrna"]:
      t[e.seqid][e.start:e.end] = e
    #efor
    return t
  #edef

  def areTandem(self, id1, id2):
    gene1 = self.genes[id1]
    gene2 = self.genes[id2]
    if gene1.seqid != gene2.seqid:
      return False
    #fi

    inregion = set([ e[-1].attr["ID"] for e in self.interval[gene1.seqid][min(gene1.end,gene2.end):max(gene1.start,gene2.start)] ])
    if len(inregion - set([gene1.attr["ID"], gene2.attr["ID"]])) == 0:
      return True
    else:
      return False
    #fi
  #edef

  def areSameStrand(self, id1, id2):
    gene1 = self.genes[id1]
    gene2 = self.genes[id2]
    return gene1.strand == gene2.strand
  #edef

###############################################################################

def readGFF3File(filename):
  G = []
  with (gzip.open(filename, "rt") if filename[-2:] == "gz" else open(filename, "r")) as gffFile:
    gffReader = csv.reader(gffFile, delimiter="\t", quotechar='"')
    for row in gffReader:
      if len(row) != 9:
        continue
      #fi
      G.append(parseGFF3entry(row))
    #efor
  #ewith
  return GFF3(G)
#edef

###############################################################################


def readMapping(mapping):
  import csv
  M = []
  with open(mapping, "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t")
    for (g1, g2) in reader:
      M.append((g1,g2))
    #efor
  #ewith
  return M
#edef
