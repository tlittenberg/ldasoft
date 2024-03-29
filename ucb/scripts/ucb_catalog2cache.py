import glob
import os

cwd = os.getcwd()
dirnames = sorted(glob.glob(cwd+"/*/catalog*/"))
print(cwd)
cache = open('ucb_catalog.cache','w')

for dir in dirnames:
    print(dir)
    entry_file=open(dir+"entries.dat","r")
    lines=entry_file.readlines()
    lines.sort()
    for line in lines:
        name,snr,evidence = line.split()
        paramfile=open(dir+name+'_params.dat')
        f0 = (paramfile.readlines()[0]).split()[0]
        print(name,f0,snr,evidence,dir,file=cache)

cache.close()
