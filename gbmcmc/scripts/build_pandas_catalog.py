# build_pandas_catalog.py
# builds a PANDAS catalog from GBMCMC output 
# modules, etc.
import argparse
import fnmatch
import os
import shutil
import glob
try: 
    import boto3
except ImportError: 
    boto3 = None
import pandas as pd
import numpy as np
import zipfile
import time
from datetime import datetime
    
# argument parsing
parser = argparse.ArgumentParser(description='Parse a set of GBMCMC catalog outputs and build a PANDAS-accessible catalog for detections and chains.')

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

parser.add_argument("-s","--spectralPattern", help = "Pattern for additional GBMCMC outputs to parse and build PANDAS DataFrame for frequency-series residuals, data, etc. If empty residual will not be processed.")

parser.add_argument("-N","--samples", help = "Number of samples argument passed to GBMCMC. Used to properly overlap adjacent segments. If empty, uses filename for AWS data or zero for local data.")






# Parse arguments
args = parser.parse_args()
v = args.verbose

bigtic = time.time()
# set up path
outDir = os.path.dirname(args.catalogFile)
catFile = os.path.basename(args.catalogFile)
tmpDir = os.path.join(outDir,'temp')
if args.aws:
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)
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
    if boto3 is None:
        raise Exception("boto3 module, required for AWS S3 link, was not imported successfully")
    if v:
        print("searching for catalog objects matching pattern %s on AWS S3 within bucket %s." % (args.pattern, args.bucket))
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

# Local directory version
else :
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
    # For AWS, first download the file and extract it
    if args.aws:
        s3 = boto3.client('s3')
        s3.download_file(args.bucket, key, 'tempcat.zip')
        
        zip = zipfile.ZipFile('tempcat.zip')
        zip.extractall() 
    # For local, change to the relevant directory
    else:
        os.chdir(key)
    
    # find the catalog directory and use the number to infer the segment
    catDir = glob.glob('catalog_*[!.sh]')
    if len(catDir)>0:
        catDir = catDir[0]

    
        if args.aws:
            # for our AWS examples, the catalog directory includes the segment
            seg = np.int(catDir[8:])
        else:
            if key[0:4]=='seg_':
                seg = np.int(key[5:])
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
                chain_cols = ['Frequency','Frequency Derivative','Amplitude','Ecliptic Longitude','coslat','cosinc','Initial Phase','Polarization','SNR','entry match', 'waveform measure']
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
                "psi" : "Initial Phase",
                "phi" : "Polarization"})

            cat_df = pd.concat([cat_df,entry_df],axis=1)
            # add the data for this segment into the list of dataframes
            cat_list.append(cat_df)

        # for AWS S3 clean up the temp directories
        if args.aws:
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
    print('Catalog complete after %0.1f/n/n' % (time.time()-bigtic))
    


# Section for building frequency-series residual data frame
rp = args.spectralPattern
if rp is None:
    if v:
        print('No spectral processing pattern specified')
else :
    medtic = time.time()
    # AWS S3 version
    if args.aws:
        if v:
            print("searching for spectral data objects matching pattern %s on AWS S3 within bucket %s." % (rp, args.bucket))
            # Build header on output file
            print("LOCAITON: AWS \r")
            print("BUCKET: %s \r" % args.bucket)
            print("PATTERN: %s \r" % rp)


        # set up S3 connection and get list of all objects in the bucket
        s3 = boto3.resource('s3')
        bucket = s3.Bucket(args.bucket)
        objs = bucket.objects.all()

        # pick out objects matching pattern
        keys = list()
        for obj in objs:
            k = obj.key
            if fnmatch.fnmatch(k,rp):
                keys.append(k)




    else :
        if v:
            print("searching local directories for spectral data matching pattern %s" % rp)
            
        keys = glob.glob(args.spectralPattern+'*')



    if v:
        print("%i spectrum directories located" % len(keys))


        
    # First set up the pandas DF that will hold all of the frequency series data
    freq_df = pd.DataFrame([],['Frequency',
                              'Median AE noise',
                              '25 AE noise', 
                              '75 AE noise',
                              '5 AE noise',
                              '95 AE noise',
                              'Median A residual',
                              '25 A residual', 
                              '75 A residual',
                              '5 A residual',
                              '95 A residual',
                              'Median E residual',
                              '25 E residual', 
                              '75 E residual',
                              '5 E residual',
                              '95 E residual',
                              'Real A data',
                              'Imag A data',
                              'Real E data',
                              'Imag E data',
                              'A residual variance',
                              'E residual variance'])

    # Write the file (which is empty now)
    freq_df.to_hdf(catOutFile, key='spectrum',mode='a')  

    # begin the loop over keys    
    for key in keys:
        tic = time.time()
        if v:
            print("Beginning processing for %s" % str(key))
        # get the number of expected samples from the file name
        if args.samples is None:
            if args.aws:
                nstr = key.split('_')
                Nsamp = np.int(nstr[len(nstr)-1].split('.')[0])/1.0
            else:
                Nsamp = 0
        else:
            Nsamp = np.int(args.samples)
            
        # for AWS, get the ZIP file and extract it to a temporary directory
        if args.aws:
            s3 = boto3.client('s3')
            s3.download_file(args.bucket, key, 'tempcat.zip')   
            zip = zipfile.ZipFile('tempcat.zip')
            zip.extractall()
        

            # find the data directory 
            rootDir = (os.path.basename(key)).split('.')[0]
            dataDir = rootDir + '/data' 
        
        # For local, we just change to the right director
        else:
            os.chdir(key)
            dataDir = 'data'


        # read the noise data 
        dfn = pd.read_csv(dataDir+'/power_noise_t0_f0.dat',
                          sep= ' ',
                          usecols=np.arange(0,6),
                          index_col = False, 
                          names=[
                              'Frequency',
                              'Median AE noise',
                              '25 AE noise', 
                              '75 AE noise',
                              '5 AE noise',
                              '95 AE noise'])


        # trim off the padding (pad length determined by length of file and specified number of samples. Keep the middle)
        N = len(dfn)
        Npad = (N-Nsamp)/2.0
        dfn.drop(N-np.arange(Npad,0,-1),inplace=True)
        dfn.drop(np.arange(0,Npad),inplace=True)


        # read the residual data 
        dfr = pd.read_csv(dataDir+'/power_residual_t0_f0.dat',
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
        Npad = (N-Nsamp)/2.0
        dfr.drop(N-np.arange(Npad,0,-1),inplace=True)
        dfr.drop(np.arange(0,Npad),inplace=True)

        # read the residual variance data 
        dfv = pd.read_csv(dataDir+'/variance_residual_t0_f0.dat',
                          sep= ' ',
                          usecols=np.arange(0,3),
                          index_col = False, 
                          names=[
                              'Frequency',
                              'A residual variance',
                              'E residual variance'])

        # trim off the padding (pad length determined by length of file and specified number of samples. Keep the middle)
        N = len(dfv)
        Npad = (N-Nsamp)/2.0
        dfv.drop(N-np.arange(Npad,0,-1),inplace=True)
        dfv.drop(np.arange(0,Npad),inplace=True)

        # read the input data 
        dfd = pd.read_csv(dataDir+'/data_0_0.dat',
                          sep= ' ',
                          usecols=np.arange(0,5),
                          index_col = False, 
                          names=[
                              'Frequency',
                              'Real A data',
                              'Imag A data',
                              'Real E data',
                              'Imag E data'])

        # trim off the padding (pad length determined by length of file and specified number of samples. Keep the middle)
        N = len(dfd)
        Npad = (N-Nsamp)/2.0
        dfd.drop(N-np.arange(Npad,0,-1),inplace=True)
        dfd.drop(np.arange(0,Npad),inplace=True)


        # concatenate the columns together to make the data frame for this segment
        df = dfn.merge(dfr)
        df = df.merge(dfv)
        df = df.merge(dfd)
        df.set_index('Frequency')

        # clean up the temp directories
        if args.aws:
            shutil.rmtree(rootDir)
            os.remove('tempcat.zip')

        # load the concatenated dataframe form disk and add on the portion for this segment
        freq_df = pd.read_hdf(catOutFile, key='spectrum',mode='r')  
        freq_df = pd.concat([freq_df,df],axis=0)

        # write the concatenated dataframe back to disk
        freq_df.to_hdf(catOutFile, key='spectrum',mode='a')

        # Print status
        if v:
            print('Completed after %0.1f' %  (time.time() - tic))

    if v:
        print('Spectral data processing completed in %0.1f\n\n' % (time.time()-medtic))
        
        
if v:
    print('Processing completed in %0.1f. Enjoy your data!\n\n'  % (time.time() - bigtic))
