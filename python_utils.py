import sys;
import shlex, subprocess;
import multiprocessing;
import os;
import pandas as pd
import numpy as np;

###############################################################################
# Running commands
# Taken from my RNA-Seq pipeline

def run_cmd(cmd, bg = False, stdin = None, stdout = None, stderr = None):
    print '[RUN CMD] Executing: %s' % cmd
    sys.stdout.flush()
    p = subprocess.Popen(shlex.split(cmd), stdin=stdin, stdout=stdout, stderr=stderr);
    if bg :
        return p
    else:
        (pid, r) = os.waitpid(p.pid, 0);
        return r;

def run_cmd_fail(cmd, stdin = None, stdout = None, stderr = None):
    print '[RUN CMD] Executing: %s' % cmd
    sys.stdout.flush()
    p = subprocess.Popen(shlex.split(cmd), stdin=stdin, stdout=stdout, stderr=stderr);
    (pid, r) = os.waitpid(p.pid, 0);
    if r != 0:
      raise RuntimeError, "Command failed: " + str(res) + "\n" + cmd;
      sys.exit(r);
    return r;

def run_shell(cmd, bg = False):
    print '[RUN SHELL] Executing: %s' % cmd
    sys.stdout.flush()
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    if bg :
        return p
    else:
        (pid, r) = os.waitpid(p.pid, 0);
        return r;

def getCommandOutput(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    res = p.communicate()[0]
    if p.returncode != 0:
        raise RuntimeError, "Command failed: " + str(res) + "\n" + cmd;
    return res

def run_par_cmds(cmd_list, max_threads=5, stdin=None, stdout=None, stderr=None, shell=False):
  
  p = [];
  i = 0;
  retval = 0;
  cmds = len(cmd_list);

  while i < cmds:
    while len(p) < max_threads and i < cmds:
      print "RUNNING: %s" % cmd_list[i]; sys.stdout.flush();
      sys.stdout.flush()
      if not shell :
        p.append( (run_cmd(cmd_list[i], bg=True, stdin=stdin, stdout=stdout, stderr=stderr),i) );
      else :
        p.append( (run_shell(cmd_list[i], bg=True),i) );
      i = i + 1;
    #ewhile

    time.sleep(0.5);

    running   = [ (j, k) for (j,k) in p if j.poll() == None ];
    completed = [ (j, k) for (j,k) in p if j.poll() != None ];

    for (j,k) in completed:
      if j.returncode != 0:
        retval = retval + j.returncode;
        print "ERROR: Failed in cmd: %s" % cmd_list[k]; sys.stdout.flush();
        sys.stdout.flush()
      else:
        print "COMPLETED: cmd : %s" % cmd_list[k]; sys.stdout.flush();
        sys.stdout.flush()
      #fi
    #efor
    p = running;
  #ewhile
  
  return retval;
#edef

def run_seq_cmds(cmd_list, stdin=None, stdout=None, stderr=None):

  for cmd in [ x for x in cmd_list if x ]:
    retval = run_cmd(cmd, stdin=stdin, stdout=stdout, stderr=stderr);
    if retval != 0:
      print "ERROR: Failed on cmd: %s" % cmd;
      sys.stdout.flush()
      return retval;
    #fi
  #efor

  return 0;
#edef

###############################################################################
  # DNA utilities

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

###############################################################################
  # Color utilities

default_colors = [(0, 0, 255), (255, 255, 255), (255, 0, 0)]  # [BLUE, WHITE, RED]

# from http://stackoverflow.com/questions/20792445/calculate-rgb-value-for-a-range-of-values-to-create-heat-map
def convert_to_rgb(minval, maxval, val, colors=default_colors):
  max_index = len(colors)-1
  v = float(val-minval) / float(maxval-minval) * max_index
  i1, i2 = int(v), min(int(v)+1, max_index)
  (r1, g1, b1), (r2, g2, b2) = colors[i1], colors[i2]
  f = v - i1
  return int(r1 + f*(r2-r1)), int(g1 + f*(g2-g1)), int(b1 + f*(b2-b1))
#edef

def rgb2frac((r,g,b)):
  return (r/255.0, g/255.0, b/255.0);
#edef

def rgb2hex((r, g, b)):
  return '%02x%02x%02x' % (r, g, b)
#edef

###############################################################################
# Natural string sorting
# From http://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside

import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def natural_sort(list):
  return sorted(list, key=natural_keys)
#edef

###############################################################################

# Scale a set of numbers in such a way that 'mid' is always in the middle. This will produce two functions
# 1: All numbers less than 'mid' are scaled to be between 0 and 0.5.
# 2: All numbers >= to 'mid' are scaled to be between 0.5 and 1.0.
def norm(min, mid, max, x):

  nmin = 0;
  nmid = 0.5
  nmax = 1.0

  if x > max:
    return nmax;
  if x < min:
    return nmin;
  if x < mid:
    return  (x-min) * (nmid / (mid - min))
  else:
    return (x-mid) * ((nmax - nmid) / (max - mid)) + nmid
  #fi
#edef

###############################################################################

def rep_to_df(rep):
  names = rep.Names;
  data = rep();

  df = pd.DataFrame({ n : d for (n,d) in zip(names, data) })

  return df[names];
#edef

