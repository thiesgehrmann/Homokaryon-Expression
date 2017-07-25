
def read_svg(svg_file):

  fd = open(svg_file, 'r');
  svg = fd.readlines();
  fd.close();

  return svg;

#edef

def write_svg(svg, out_file):

  fd = open(out_file, 'w');
  [ fd.write(x) for x in svg ]
  fd.close();

#edef

def replace(svg, orig, new):

  return [ l.replace(orig, new) for l in svg ]

#edef

def replacem(svg, d):

  for o in d.keys():
    svg = replace(svg, o, d[o]);
  #efor
  return svg;

#edef
