# build_pandas_catalog.py
# builds a PANDAS catalog from GBMCMC output 
# modules, etc.
import argparse
import fnmatch
import os
import shutil
import glob
import boto3
import pandas as pd
import numpy as np
import zipfile
import time
from datetime import datetime
    
# argument parsing
parser = argparse.ArgumentParser(description='Build a list of GBMCMC catalog objects that can be used by make_panda_cat.py to build a PANDAS-accessible catalog.')

parser.add_argument("pattern", help="the pattern to search for in either the S3 buckets or the local directory")

parser.add_argument("catalogFile", help = "file path for the primary catalog output. Auxilliary files will be placed in the same directory")


parser.add_argument("-v", "--verbose", help="provide additional information at the console output",
                    action="store_true")


parser.add_argument("-a", "--aws", help="search for the objects in an AWS S3 bucket. If false, pattern is assumed to point to a local directory",
                    action="store_true")

parser.add_argument("-b", "--bucket", help="for AWS S3 data, name of the S3 bucket in which to search for pattern.")

parser.add_argument("-n", "--name", help="Name for the catalog (different than the filename)")

parser.add_argument("-p","--parent", help = "Name for the parent catalog (if it exists. For use in tracking sources across catalogs)")

parser.add_argument("-d","--duration", help = "Observation duration of catalog (in seconds)")




# Parse arguments
args = parser.parse_args()
v = args.verbose

bigtic = time.time()
# set up path
outDir = os.path.dirname(args.catalogFile)
catFile = os.path.basename(args.catalogFile)
tmpDir = os.path.join(outDir,'temp')
os.chdir(tmpDir)

# Build metadata 
catName = args.name # name of the catalog
metaDict = {};  # initialize the dictionary that will be used for metadata (will be converted to Pandas DF). You can add whatever you want.
metaDict['Observation Time'] = np.float(args.duration) #observation time
metaDict['parent'] = args.parent # name of the parent catalog (if it exists)
metaDict['Build Time'] = datetime.strftime(datetime.now(),'%Y-%m-%d %H:%M:%S') # creation date for catalog
if v:
    print("Preparing catalog with the following metadata:\r")
    print("NAME: %s\r" % catName)
    print("OBSERVATION TIME: %s sec\r" % metaDict['Observation Time'])
    print("PARENT CATALOG: %s\r" % metaDict['parent'])
    print("BUILD TIME: %s" % metaDict['Build Time'])


# AWS S3 version
if args.aws:
    if v:
        print("searching for objects matching pattern %s on AWS S3 within bucket %s." % (args.pattern, args.bucket))
        # Build header on output file
        print("LOCAITON: AWS \r")
        print("BUCKET: %s \r" % args.bucket)
        print("PATTERN: %s \r" % args.pattern)
       

    # set up S3 connection and get list of all objects in the bucket
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(args.bucket)
    objs = bucket.objects.all()
    
    # pick out objects matching pattern
    pat = args.pattern
    keys = list()
    for obj in objs:
        k = obj.key
        if fnmatch.fnmatch(k,pat):
            keys.append(k)


      
    
else :
    print("searching local directories for objects matching pattern %s" % args.pattern)
    


if v:
    print("%i segment catalogs located" % len(keys))

    
# build the metadata data Frame and write to the file
metaDF = pd.DataFrame(data = metaDict, index={catName})
catOutFile = os.path.join(outDir,catFile)
metaDF.to_hdf(catOutFile, key='metadata',mode='w')  # this overwrites so you get a new catalog if you run this more than once

# Iterate over all the keys to make catalog (takes a little while)
cat_list = list()
for key in keys:
    tic = time.time()
    if v:
        print("Beginning processing for %s" % key)
    # first get the file and extract it
    if args.aws:
        s3 = boto3.client('s3')
        s3.download_file(args.bucket, key, 'tempcat.zip')
        
    zip = zipfile.ZipFile('tempcat.zip')
    zip.extractall()  
    
    # find the catalog directory and use the number to infer the segment
    catDir = glob.glob('catalog_*')[0]  # what if there is more than 1?
    seg = np.int(catDir[8:])

    # check if there is anything in the entries file
    catSize = os.path.getsize(catDir+'/'+'entries.dat')
    if catSize > 1:
        # read the entries file
        cat_df = pd.read_table(catDir+'/'+'entries.dat',delimiter = ' ',index_col=0,names=['SNR','evidence'])
        cat_df['segment']=seg
        segFile = "%s_chains_%is.h5" % (os.path.splitext(catFile)[0],100*np.floor(seg/100))
        segOutFile = os.path.join(outDir,segFile)
        cat_df['chain file']=segFile
        
        # read the history file (if present) and concatenate it with the entries file
        histFile = catDir+'/'+'history.dat'
        if os.path.isfile(histFile):
            if  (os.path.getsize(histFile) > 2 ):
                hist_df = pd.read_table(histFile,delimiter = ' ', index_col=1,names=['parent'])
                hist_df = hist_df[~hist_df.index.duplicated(keep='first')]
                cat_df = pd.concat([cat_df,hist_df],axis=1)
            else :
                cat_df['parent']=''
        else:
            cat_df['parent']=''

        # read each entry to get the median parameters for the point-estimate data frame & store the chains for the other data frames
        entryNameList = list(cat_df.index)
        dfCols = ['f','fdot','amp','lon','coslat','cosinc','psi','phi'] # this is hardcoded, if the chain files had headers, we could make this smart
        dfs = list()
        for ii,entryName in enumerate(entryNameList):
            df = pd.read_table(catDir + '/' + entryName +'_params.dat',delimiter = ' ',index_col = False, names=dfCols)
            df['name']=entryName
            df=df.set_index('name')
            dfs.append(df)

            # Build data frame for chain for this particular source
            chain_cols = ['Frequency','Frequency Derivative','Amplitude','Ecliptic Longitude','coslat','cosinc','Initial Phase','Polarization','SNR','entry match', 'waveform measure'];
            chain_df = pd.read_table(catDir+'/'+entryName +'_chain.dat',delimiter = ' ',index_col = False, names=chain_cols)
            chain_df['Ecliptic Latitude']=np.arccos(chain_df['coslat'])-np.pi/2.0  # I think the conventions are correct, should check
            chain_df['Inclination']=np.arccos(chain_df['cosinc']) - np.pi/2.0 # I think the conventions are correct, should check
            #chain_df = chain_df.drop(columns=['coslat','cosinc'])

            # write chain for this entry to the HDF5 file
            
            chain_df.to_hdf(segOutFile, key=entryName + '_chain',mode='a')

        # combine the parameters and entries catalogs
        entry_df = pd.concat(dfs)
        entry_df['Ecliptic Latitude']=np.arccos(entry_df['coslat'])-np.pi/2.0
        entry_df['Inclination']=np.arccos(entry_df['cosinc']) - np.pi/2.0
        #entry_df = entry_df.drop(columns=['coslat','cosinc'])
        entry_df = entry_df.rename(columns={
            "f" : "Frequency", 
            "fdot" : "Frequency Derivative",
            "amp" : "Amplitude",
            "lon" : "Ecliptic Longitude",
            "psi" : "Initial Phase",
            "phi" : "Polarization"})

        cat_df = pd.concat([cat_df,entry_df],axis=1)
        # add the data for this segment into the list of dataframes
        cat_list.append(cat_df)

    # clean up the temp directories
    shutil.rmtree(catDir)
    os.remove('tempcat.zip')
    
    # Print status
    if v:
        print('Completed after %0.1f' %  (time.time() - tic))


        

# concatenate the data frames for all segments into one data frame and write it to the HDF5 file
totcat_df = pd.concat(cat_list)
totcat_df = totcat_df.replace(np.nan, '', regex=True)
totcat_df.to_hdf(catOutFile, key='detections',mode='a')

if v:
    print('Catalog processing completed in %0.1f '  % (time.time() - bigtic))