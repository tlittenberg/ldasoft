import glob
import os
import subprocess
import numpy as np

rundir=os.getcwd()
ucb_dirs = sorted(glob.glob("ucb/seg*"))
detections=[]
frequency=[]

#Post process UCB model
ps=[]
for dir in ucb_dirs:
    os.chdir(dir)
    
    fbins = np.loadtxt('data/power_data_0.dat',usecols=(0))
    med_f = np.median(fbins)
    
    frequency.append(med_f)


    evidence = np.loadtxt('evidence.dat',usecols=(1))
    max_N = np.argmax(evidence)
    
    detections.append(max_N)
    
    #os.system('bash example_gb_catalog.sh ' + str(max_N))
    p=subprocess.Popen('bash example_gb_catalog.sh ' + str(max_N) + ' > post.out', shell=True)
    ps.append(p)  
    os.chdir(rundir)

for p in ps:
    p.wait()

#Save UCBs to cache file
dirnames = glob.glob("ucb/seg*/catalog*/")
cache = open(rundir+'/gb_catalog.cache','w')

for dir in dirnames:
    entry_file=open(dir+"entries.dat","r")
    lines=entry_file.readlines()
    for line in lines:
        name,snr,evidence = line.split()
        paramfile=open(dir+name+'_params.dat')
        f0 = (paramfile.readlines()[0]).split()[0]
        print(name,f0,snr,evidence,dir,file=cache)

cache.close()


#Post process VGBs
vgb_dirs = sorted(glob.glob(rundir+"/vgb/seg*"))
detections=[]
frequency=[]

ps=[]
for dir in vgb_dirs:
    os.chdir(dir)
    fbins = np.loadtxt('data/power_data_0.dat',usecols=(0))
    med_f = np.median(fbins)
    min_f = np.min(fbins)
    
    frequency.append(med_f)

    max_N = 1
    
    detections.append(max_N)
    
    p=subprocess.Popen('bash ../example_gb_catalog.sh ' + str(max_N) + ' ' + str(min_f) + ' > post.out', shell=True)
    ps.append(p)
    os.chdir(rundir)

for p in ps:
    p.wait()
