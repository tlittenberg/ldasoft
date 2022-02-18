# build_pandas_catalog.py
# builds a PANDAS catalog from MBH output 
# modules, etc.
import argparse
import fnmatch
import os
import shutil
import glob
import pandas as pd
import numpy as np
import time
from datetime import datetime
    
# argument parsing
parser = argparse.ArgumentParser(description='Parse a set of MBH catalog outputs and build a PANDAS-accessible catalog for detections and chains.')

parser.add_argument("pattern", help="Search pattern for the directories containing the sources")

parser.add_argument("catalogFile", help = "file path for the primary catalog output. Auxilliary files will be placed in the same directory")

parser.add_argument("-v", "--verbose", help="provide additional information at the console output",
                    action="store_true")

parser.add_argument("-n", "--name", help="Name for the catalog (different than the filename)")

parser.add_argument("-p","--parent", help = "Name for the parent catalog (if it exists. For use in tracking sources across catalogs)")

parser.add_argument("-d","--duration", help = "Observation duration of catalog (in seconds)")


# Parse arguments
args = parser.parse_args()
v = args.verbose
dur = np.float(args.duration)

bigtic = time.time()
# set up path
outputFile = args.catalogFile
outDir = os.path.dirname(outputFile)
catFile = os.path.basename(outputFile)

# Build metadata 
catName = args.name # name of the catalog
metaDict = {};  # initialize the dictionary that will be used for metadata (will be converted to Pandas DF). You can add whatever you want.
metaDict['Observation Time'] = dur #observation time
metaDict['parent'] = args.parent # name of the parent catalog (if it exists)
metaDict['Build Time'] = datetime.strftime(datetime.now(),'%Y-%m-%d %H:%M:%S') # creation date for catalog
if v:
    print("Preparing catalog with the following metadata:\r")
    print("NAME: %s\r" % catName)
    print("OBSERVATION TIME: %s sec\r" % metaDict['Observation Time'])
    print("PARENT CATALOG: %s\r" % metaDict['parent'])
    print("BUILD TIME: %s" % metaDict['Build Time'])


if v:
	print("searching local directories for objects matching pattern %s" % args.pattern)

keys = glob.glob(args.pattern+'*')
    
if v:
    print("%i sources located" % len(keys))

    
# build the metadata data Frame and write to the file
metaDF = pd.DataFrame(data = metaDict, index={catName})
catOutFile = os.path.join(outDir,catFile)
metaDF.to_hdf(catOutFile, key='metadata',mode='w')  # this overwrites so you get a new catalog if you run this more than once

colNames = ['Iteration', 
            'Log Likelihood', 
            'Mass 1', 
            'Mass 2', 
            'Spin 1',
            'Spin 2', 
            'Merger Phase', 
            'Barycenter Merge Time', 
            'Luminosity Distance', 
            'cos ecliptic colatitude', 
            'Ecliptic Longitude',
            'Polarization',
            'cos inclination',
            'Index 2',
            'Theta L',
            'phi_L',
            'A noise',
            'E noise'];

# Iterate over all the keys to make catalog (takes a little while)
cat_list = list()
pointEstimates = np.empty([len(keys),len(colNames)-1])
srcs = list()

for key in keys:
    tic = time.time()
    if v:
        print("Reading source %s" % os.path.basename(key))
    
    # read chain file
    chainDF = pd.read_csv(key + '/chain.dat', sep = ' ', index_col=0,names = colNames)
    chainDF['Ecliptic Latitude'] = np.pi/2-np.arccos(np.array(chainDF['cos ecliptic colatitude']))

    # take the median to get the point-estimates (maybe not good for multi-modal?)
    src = chainDF.median()
    src['name'] = 'MBH%09i' % src['Barycenter Merge Time']
    src['chain file'] = os.path.basename(outputFile)
    src['source number'] = int(os.path.basename(key)[-4:])
    src['Observation Time'] = dur
    srcs.append(src)

    chainDF.to_hdf(outputFile,key=('%s_chain' % src['name']), mode = 'a')
    
    # read the power reconstruction data 
    powerDF = pd.read_csv(key + '/power_reconstruction.dat',
        sep= ' ',
        usecols=np.arange(0,11),
        index_col = False, 
        names=[
            'Frequency',
            'Median A power',
            '25 A power', 
            '75 A power',
            '5 A power',
            '95 A power',
            'Median E power',
            '25 E power', 
            '75 E power',
            '5 E power',
            '95 E power'])

    powerDF.to_hdf(outputFile,key=('%s_power' % src['name']), mode = 'a')
    
    
catDF = pd.DataFrame(srcs)
catDF.set_index('name', inplace=True)
catDF.to_hdf(outputFile,key='detections',mode='a')

if v:
    print('Processing completed in %0.1f. Enjoy your data!\n\n'  % (time.time() - bigtic))
