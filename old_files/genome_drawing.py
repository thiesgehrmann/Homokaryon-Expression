#!/usr/bin/env python

"""
Modified from https://gist.githubusercontent.com/daler/c98fc410282d7570efc3/raw/ea7f2736909ae73958c14c5be7f31ae7c65a9348/ideograms.py
Demonstrates plotting chromosome ideograms and genes (or any features, really)
using matplotlib.

1) Assumes a file from UCSC's Table Browser from the "cytoBandIdeo" table,
saved as "ideogram.txt". Lines look like this::

    #chrom  chromStart  chromEnd  name    gieStain
    chr1    0           2300000   p36.33  gneg
    chr1    2300000     5300000   p36.32  gpos25
    chr1    5300000     7100000   p36.31  gneg

2) Assumes another file, "ucsc_genes.txt", which is a BED format file
   downloaded from UCSC's Table Browser. This script will work with any
   BED-format file.

"""

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas


# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
#
def chromosome_collections(df, y_positions, height,  **kwargs):
    """

    Yields BrokenBarHCollection of features that can be added to an Axes
    object.

    Parameters
    ----------

    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.

    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection

    height : float
        Height of each BrokenBarHCollection

    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print chrom
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']
#edef

def draw_genome(chromosome_list, ideogram, ideogene, outfig, ideogroup=None):

  # Height of each ideogram
  chrom_height = 1
  
  # Spacing between consecutive ideograms
  chrom_spacing = 0.6
  
  # Height of the gene track. Should be smaller than `chrom_spacing` in order to
  # fit correctly
  gene_height  = 1.0
  group_height = 0.3
  
  # Padding between the top of a gene track and its corresponding ideogram
  gene_padding  = 0.1
  group_padding = 0.1
  
  # Width, height (in inches)
  figsize = (10,5)
  
  # Decide which chromosomes to use
  chromosome_list = chromosome_list
  
  # Keep track of the y positions for ideograms and genes for each chromosome,
  # and the center of each ideogram (which is where we'll put the ytick labels)
  ybase = 0
  chrom_ybase = {}
  gene_ybase = {}
  group_ybase = {}
  chrom_centers = {}
  
  # Iterate in reverse so that items in the beginning of `chromosome_list` will
  # appear at the top of the plot
  for chrom in chromosome_list[::-1]:
      chrom_ybase[chrom] = ybase
      chrom_centers[chrom] = ybase + chrom_height / 2.
      gene_ybase[chrom] = ybase # - gene_height - gene_padding
      group_ybase[chrom] = ybase - group_height - group_padding
      ybase += chrom_height + chrom_spacing
  
  # Read in ideogram.txt, downloaded from UCSC Table Browser
  ideo = ideogram;
  ideo.columns = ['chrom', 'start', 'end', 'name', 'colors']
  
  # Filter out chromosomes not in our list
  ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]
  
  # Add a new column for width
  ideo['width'] = ideo.end - ideo.start
  
  # Same thing for genes
  genes = ideogene
  genes.columns = ['chrom', 'start', 'end', 'name', 'mean', 'colors']
  genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
  genes['width'] = genes.end - genes.start

  if ideogroup is not None:
    groups = ideogroup
    groups.columns = ['chrom', 'start', 'end', 'colors']
    groups = groups[groups.chrom.apply(lambda x: x in chromosome_list)]
    groups['width'] = genes.end - genes.start

  fig = plt.figure(figsize=figsize)
  ax = fig.add_subplot(111)
  
  # Now all we have to do is call our function for the ideogram data...
  print("adding ideograms...")
  for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
      ax.add_collection(collection)
  
  # ...and the gene data
  print("adding genes...")
  for collection in chromosome_collections(genes, gene_ybase, gene_height, alpha=0.5, linewidths=0.05):
      ax.add_collection(collection)

  if ideogroup is not None:
    print ideogroup
    print("adding groups...")                
    for collection in chromosome_collections(groups, group_ybase, group_height, alpha=0.5, linewidths=0):
      ax.add_collection(collection)

  # Axes tweaking
  ax.set_yticks([chrom_centers[i] for i in chromosome_list])
  ax.set_yticklabels(chromosome_list)
  ax.axis('tight')
  plt.savefig(outfig);
#edef

