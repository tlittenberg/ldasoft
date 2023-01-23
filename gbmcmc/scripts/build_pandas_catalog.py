# build_pandas_catalog.py
# builds a PANDAS catalog from GBMCMC output 
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
parser = argparse.ArgumentParser(description='Parse a set of GBMCMC catalog outputs and build a PANDAS-accessible catalog for detections and chains.')

parser.add_argument("pattern", help="the pattern to search for UCB segment output directories")

parser.add_argument("catalogFile", help = "file path for the primary catalog output. Auxilliary files will be placed in the same directory")

parser.add_argument("-v", "--verbose", help="provide additional information at the console output",
                    action="store_true")

parser.add_argument("-n", "--name", help="Name for the catalog (different than the filename)")

parser.add_argument("-p","--parent", help = "Name for the parent catalog (if it exists. For use in tracking sources across catalogs)")

parser.add_argument("-d","--duration", help = "Observation duration of catalog (in seconds)")

parser.add_argument("-N","--samples",default = 32,  help = "Number of samples argument passed to GBMCMC. Used to properly overlap adjacent segments. Default = 32")


# Parse arguments
args = parser.parse_args()
v = args.verbose

bigtic = time.time()
# set up path
outDir = os.path.dirname(args.catalogFile)
catFile = os.path.basename(args.catalogFile)
tmpDir = os.path.join(outDir,'temp')

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

# Local directory version

if v:
    print("searching local directories for objects matching pattern %s" % args.pattern)
    
keys = glob.glob(args.pattern+'*')
    

if v:
    print("%i catalogs located" % len(keys))

    
# build the metadata data Frame and write to the file
metaDF = pd.DataFrame(data = metaDict, index={catName})
catOutFile = os.path.join(outDir,catFile)
metaDF.to_hdf(catOutFile, key='metadata',mode='w')  # this overwrites so you get a new catalog if you run this more than once

# Iterate over all the keys to make catalog (takes a little while)
cat_list = list()
ss = 0;
for key in keys:
    tic = time.time()
    if v:
        print("Beginning processing for %s" % key)

    os.chdir(key)
    
    # find the catalog directory and use the number to infer the segment
    catDir = glob.glob('catalog_*[!.sh]')
    if len(catDir)>0:
        catDir = catDir[0]

        if key[0:4]=='seg_':
            seg = np.int(key[-3:])
        else:
            seg = ss;
            ss = ss+1;

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
                    cat_df = cat_df.join(hist_df, how='left')
                    cat_df.replace(np.nan,'',regex=True,inplace=True)
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
                chain_cols = ['Frequency','Frequency Derivative','Amplitude','Ecliptic Longitude','coslat','cosinc','Polarization','Initial Phase','SNR','entry match', 'waveform measure']
                chain_df = pd.read_table(catDir+'/'+entryName +'_chain.dat',delimiter = ' ',index_col = False, names=chain_cols)
                chain_df['Ecliptic Latitude']=np.pi/2.0-np.arccos(pd.to_numeric(chain_df['coslat'],'coerce'))
                chain_df['Inclination']=np.arccos(pd.to_numeric(chain_df['cosinc'],'coerce'))
                #chain_df = chain_df.drop(columns=['coslat','cosinc'])

                # write chain for this entry to the HDF5 file
    
                chain_df.to_hdf(segOutFile, key=entryName + '_chain',mode='a')
    
                # Build data frame for waveform for this particular source
                wave_cols = ['Frequency','Waveform Real A','Waveform Imag A','Waveform Real E','Waveform Imag E']
                wave_df = pd.read_table(catDir+'/'+entryName +'_waveform.dat',delimiter = ' ',index_col = False, names=wave_cols)
                # write chain for this entry to the HDF5 file
    
                wave_df.to_hdf(segOutFile, key=entryName + '_wave',mode='a')

            # combine the parameters and entries catalogs
            entry_df = pd.concat(dfs)
            entry_df['Ecliptic Latitude']=np.pi/2.0-np.arccos(pd.to_numeric(entry_df['coslat'],'coerce'))
            entry_df['Inclination']=np.arccos(pd.to_numeric(entry_df['cosinc'],'coerce'))
            #entry_df = entry_df.drop(columns=['coslat','cosinc'])
            entry_df = entry_df.rename(columns={
                "f" : "Frequency", 
                "fdot" : "Frequency Derivative",
                "amp" : "Amplitude",
                "lon" : "Ecliptic Longitude",
                "psi" : "Polarization",
                "phi" : "Initial Phase"})

            cat_df = pd.concat([cat_df,entry_df],axis=1)
            # add the data for this segment into the list of dataframes
            cat_list.append(cat_df)
            
        # Print status
        if v:
            print('Completed after %0.1f s' %  (time.time() - tic))


        

# concatenate the data frames for all segments into one data frame and write it to the HDF5 file
totcat_df = pd.concat(cat_list)
totcat_df = totcat_df.replace(np.nan, '', regex=True)
totcat_df.to_hdf(catOutFile, key='detections',mode='a')


if v:
    print('Catalog complete after %0.1f s\n\n' % (time.time()-bigtic))
    


# Section for building frequency-series residual data frame
medtic = time.time()
if v:
    print("searching local directories for spectral data matching pattern %s" % args.pattern)

keys = glob.glob(args.pattern+'*')

if v:
    print("%i spectrum directories located" % len(keys))

# begin the loop over keys    
for key in keys:
    tic = time.time()
    if v:
        print("Beginning processing for %s" % str(key))
    
    os.chdir(key)
    dataDir = 'data'

    # read the residual data 
    dfr = pd.read_csv(dataDir+'/power_residual_t0.dat',
        sep= ' ',
        usecols=np.arange(0,11),
        index_col = False, 
        names=[
            'Frequency',
            'Median A residual',
            '25 A residual', 
            '75 A residual',
            '5 A residual',
            '95 A residual',
            'Median E residual',
            '25 E residual', 
            '75 E residual',
            '5 E residual',
            '95 E residual'])

    # trim off the padding (pad length determined by length of file and specified number of samples. Keep the middle)
    N = len(dfr)
    Npad = np.int(args.samples)
    dfr.drop(N-np.arange(Npad,0,-1),inplace=True)
    dfr.drop(np.arange(0,Npad),inplace=True)

    # read the residual variance data 
    dfv = pd.read_csv(dataDir+'/variance_residual_t0.dat',
        sep= ' ',
        usecols=np.arange(0,3),
        index_col = False, 
        names=[
            'Frequency',
            'A residual variance',
            'E residual variance'])

    # trim off the padding (pad length determined by length of file and specified number of samples. Keep the middle)
    Npad = np.int(args.samples)
    dfv.drop(N-np.arange(Npad,0,-1),inplace=True)
    dfv.drop(np.arange(0,Npad),inplace=True)
    
    # read the power reconstruction data 
    dfp = pd.read_csv(dataDir+'/power_reconstruction_t0.dat',
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

    # trim off the padding (pad length determined by length of file and specified number of samples. Keep the middle)
    dfp.drop(N-np.arange(Npad,0,-1),inplace=True)
    dfp.drop(np.arange(0,Npad),inplace=True)
     
     
    # concatenate the columns together to make the data frame for this segment
    df = dfr.merge(dfv)
    df = df.merge(dfp)
    df.set_index('Frequency')


    # load the concatenated dataframe form disk and add on the portion for this segment
    try:
        freq_df = pd.read_hdf(catOutFile, key='spectrum',mode='r')  
        freq_df = pd.concat([freq_df,df],axis=0)
    except:
        freq_df = df
            
    # write the concatenated dataframe back to disk
    freq_df.to_hdf(catOutFile, key='spectrum',mode='a')

    # Print status
    if v:
        print('Completed after %0.1f s' %  (time.time() - tic))

if v:
    print('Spectral data processing completed in %0.1f s \n\n' % (time.time()-medtic))

if v:
    print('Processing completed in %0.1f s. Enjoy your data! \n\n'  % (time.time() - bigtic))
