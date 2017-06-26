#!/usr/bin/env python

import numpy as np
from matplotlib import rcParams
import corner
import sys, getopt

def main(argv):
  inputfile = ''
  outputfile = ''
  try:
    opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
  except getopt.GetoptError:
    print 'test.py -i <inputfile> -o <outputfile>'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print 'test.py -i <inputfile> -o <outputfile>'
      sys.exit()
    elif opt in ("-i", "--ifile"):
      inputfile = arg
    elif opt in ("-o", "--ofile"):
      outputfile = arg
  print 'Input file is "', inputfile
  print 'Output file is "', outputfile


  data = np.loadtxt(inputfile, usecols=(0,1,2))

  rcParams.update({'font.size': 8})

  figure = corner.corner(data, quantiles=[0.1, 0.5, 0.9], show_titles=True, title_kwargs={"fontsize": 8},labels=[r"$f$ [mHz]", r"$\log_{10}\frac{df}{dt}$ [s$^{-2}$]",r"$\log_{10}\frac{d^2f}{dt^2}$ [s$^{-3}$]"])

  figure.savefig(outputfile, dpi=300)


if __name__ == "__main__":
  main(sys.argv[1:])
